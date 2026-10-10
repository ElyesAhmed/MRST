% optimizeModel2D - optimize NPV for example-model of this folder

mrstModule add ad-core ad-blackoil ad-props optimization spe10
setupModel2D

% Create model-object of class TwoPhaseOilWaterModel
model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
% Set initial state and run simulation:
state0 = initResSol(G, 200*barsa, [0, 1]);

problem     = packSimulationProblem(state0, model, schedule, 'tmp', 'Name', 'model2D', 'Description', 'model2d');

numRealizations = 1;
problems = cell(1, numRealizations);
for k = 1:numRealizations
    problems{k} = problem;
    problems{k}.seed = k;
end
samples = struct('problem', {problems}, 'num', numRealizations);
%% Set up box limits for scaling and define function evaluation handle
G = model.G;

% Setup target and limit bounds
% Initialize bounds structure
bounds = [];

for j = 1:length(schedule.control)    
    W = schedule.control(j).W;    
    for i = 1:length(W)
        if W(i).sign > 0  % Injector
            % Get reference rate (use first control step as reference if needed)
            ref_rate = W(i).val;            
            % Define bounds for injector
            step_bounds = [0, ref_rate*3];  % Rate bounds   
            bounds = [bounds;step_bounds];
        end
    end


    for i = 1:length(W)
        if W(i).sign < 0  % Producer
            % Define bounds for producer

            ref_bhp = W(i).val/barsa;            
            step_bounds = [100, 250] * barsa;                % BHP bounds
            bounds = [bounds;step_bounds];

        end
    end
end


bnds = processBounds(W, ...
    'rate(inj)', bounds(1,:), ...     % allowed target rates
    'rate(inj)', bounds(2,:), ...     % allowed target rates
    'bhp(prod)', bounds(3,:), ...     % allowed target rates
    'bhp(prod)', bounds(4,:),  ...
    'rate(inj)', bounds(5,:), ...     % allowed target rates
    'rate(inj)', bounds(6,:), ...     % allowed target rates
    'bhp(prod)', bounds(7,:), ...     % allowed target rates
    'bhp(prod)', bounds(8,:), ...
    'rate(inj)', bounds(9,:), ...     % allowed target rates
    'rate(inj)', bounds(10,:), ...     % allowed target rates
    'bhp(prod)', bounds(11,:), ...     % allowed target rates
    'bhp(prod)', bounds(12,:), ...
    'rate(inj)', bounds(13,:), ...     % allowed target rates
    'rate(inj)', bounds(14,:), ...     % allowed target rates
    'bhp(prod)', bounds(15,:), ...     % allowed target rates
    'bhp(prod)', bounds(16,:));     % allowed target bhp

[maps, u_init] = setupSimulationControlMappings(schedule, bnds);

[maps, u_init] = reduceControlMappings(maps, u_init, 'exludeLims', true);
u_init = [    0.8854
    0.6368
    1.0000
         0
    0.4819
    0.6528
    1.0000
         0
    0.4598
    0.6173
    0.9484
    0.0739
    0.6386
    0.6028
    0.8171
    0.2010];

objectiveScaling = 3.2e7;
% Define objective function
d   = 0.0;    % yearly discount factor
ro  = 60;      % oil revenue/price ($/stb)
rwp =  6;      % water production handling costs ($/stb)
rwi =  6;      % water injection cost ($/stb) 
npvopts = {'OilPrice',             ro , ...
           'WaterProductionCost', rwp , ...
           'WaterInjectionCost',  rwi , ...
           'DiscountFactor',        d};


%% Define objective functions
objectiveOptsDiscounted = {'OilPrice', ro, ...
                          'WaterInjectionCost', rwi, ...
                          'WaterProductionCost', rwp, ...
                          'DiscountFactor', 0.5};

objectiveOptsUndiscounted = {'OilPrice', ro, ...
                            'WaterInjectionCost', rwi, ...
                            'WaterProductionCost', rwp, ...
                            'DiscountFactor', 0};

objFunDiscounted = @(model, states, schedule, varargin) ...
    NPVOW(model, states, schedule, varargin{:}, objectiveOptsDiscounted{:});

objFunUndiscounted = @(model, states, schedule, varargin) ...
    NPVOW(model, states, schedule, varargin{:}, objectiveOptsUndiscounted{:});

%% Option 1: Nash-based Alternating Optimization
fprintf('\n=== Option 1: Nash-based Alternating Optimization ===\n');

% Define control masks for each objective
% Assuming controls order: [bhp(inj), rate(inj), bhp(prod), lrat(prod)]
inj_controls_mask = [false, true];   % Injection well controls: [bhp(inj), rate(inj)]
prod_controls_mask = [true, false];  % Production well controls: [bhp(prod), lrat(prod)]

% Setup separate mappings for each objective
%[maps_disc, u_disc_init] = setupSimulationControlMappings(schedule, bnds, inj_controls_mask);
%[maps_undisc, u_undisc_init] = setupSimulationControlMappings(schedule, bnds, prod_controls_mask);
[maps_disc, u_disc_init] = reduceControlMappings(maps, u_init, 'type', {'bhp', 'lrat'});
[maps_undisc, u_undisc_init] = reduceControlMappings(maps, u_init, 'type', {'rate', 'rate'});

% Setup optimization problems
p_discounted = OptimizationProblem(samples, ...
    'name', 'Nash_Discounted', ...
    'objective', struct('function', objFunDiscounted, 'scaling', objectiveScaling), ...
    'maps', maps_disc, ...
    'setupType', 'simulation', ...
    'verboseSimulation', true);

p_undiscounted = OptimizationProblem(samples, ...
    'name', 'Nash_Undiscounted', ...
    'objective', struct('function', objFunUndiscounted, 'scaling', objectiveScaling), ...
    'maps', maps_undisc, ...
    'setupType', 'simulation', ...
    'verboseSimulation', true);

% Initialize controls
u_disc_current = u_disc_init;
u_undisc_current = u_undisc_init;

% Nash optimization parameters
maxNashIterations = 10;
convergence_tolerance = 1e-4;
nashResults = struct();

% Initialize problem for Nash iterations
problem_current = problems{1};
u_cc_undisc =  p_undiscounted.unscaleVariables([u_init(1:2);u_init(5:6);u_init(9:10);u_init(13:14)]);
u_cc_disc =  p_discounted.unscaleVariables([u_init(3:4);u_init(7:8);u_init(11:12);u_init(15:16)]);
problem_current_cc = p_discounted.updateProblemFun(problem_current,u_cc_disc);
problem_current_cc = p_undiscounted.updateProblemFun(problem_current_cc,u_cc_undisc);
problem_current = p_undiscounted.updateProblemFun(problem_current_cc,u_cc_undisc);
for nashIter = 1:maxNashIterations
    fprintf('\n--- Nash Iteration %d/%d ---\n', nashIter, maxNashIterations);
    
    % Phase 1: Optimize injection controls for discounted NPV
    fprintf('Optimizing injection controls (discounted NPV)...\n');
    p_discounted.reset('prompt', false)
    [u_disc_opt, history_disc] = p_discounted.optimize(...
        problem_current, 'stepInit', 0.1, 'maxIt', 2, 'regularizationParam', 0.025);
    
    % Update problem with optimized injection controls

    u_disc_current =  p_discounted.unscaleVariables(history_disc.u{1});
    problem_current = p_discounted.updateProblemFun(problem_current, 0.5*u_disc_opt+0.5*u_disc_current);    
    % Phase 2: Optimize production controls for undiscounted NPV
    fprintf('Optimizing production controls (undiscounted NPV)...\n');
    p_undiscounted.reset('prompt', false)
    [u_undisc_opt, history_undisc] = p_undiscounted.optimize(...
        problem_current, 'stepInit', 0.1, 'maxIt', 20, 'regularizationParam', 0.025);
    
    % Update problem with optimized production controls

    u_undisc_current =  p_undiscounted.unscaleVariables(history_undisc.u{1});

    problem_current = p_undiscounted.updateProblemFun(problem_current, 0.5*u_undisc_opt+0.5*u_undisc_current);    
    % Simulate with current controls
    model = problem_current.SimulatorSetup.model;
    schedule_current = problem_current.SimulatorSetup.schedule;
    [~, states] = simulateScheduleAD(state0, model, schedule_current);
    
    % Calculate both NPVs
    npv_disc = objFunDiscounted(model, states, schedule_current);
    npv_undisc = objFunUndiscounted(model, states, schedule_current);
    
    % Store results
    nashResults.iteration(nashIter) = nashIter;
    nashResults.npv_discounted(nashIter) = sum([npv_disc{:}]);
    nashResults.npv_undiscounted(nashIter) = sum([npv_undisc{:}]);
    nashResults.combined_npv(nashIter) = 0.5 * (nashResults.npv_discounted(nashIter) + nashResults.npv_undiscounted(nashIter));
    
    fprintf('NPV - Discounted: $%.2f, Undiscounted: $%.2f, Combined: $%.2f\n', ...
            nashResults.npv_discounted(nashIter), ...
            nashResults.npv_undiscounted(nashIter), ...
            nashResults.combined_npv(nashIter));
    
%     % Check convergence
%     if nashIter > 1
%         improvement = abs(nashResults.combined_npv(nashIter) - nashResults.combined_npv(nashIter-1)) / ...
%                      abs(nashResults.combined_npv(nashIter-1));
%         fprintf('Improvement: %.6f\n', improvement);
%         if improvement < convergence_tolerance
%             fprintf('Convergence achieved after %d iterations\n', nashIter);
%             break;
%         end
%     end
end

%% Helper Functions (Place these at the end of your file)
function full_controls = combineNashControls(u_disc, u_undisc, maps_disc, maps_undisc, maps_full)
    % Convert to physical space
    u_disc_phys = maps_disc.unscaleVariables(u_disc);
    u_undisc_phys = maps_undisc.unscaleVariables(u_undisc);
    
    % Combine based on control structure
    full_controls_phys = zeros(1, 4);
    full_controls_phys(1:2) = u_disc_phys;    % Injection controls
    full_controls_phys(3:4) = u_undisc_phys;  % Production controls
    
    % Convert back to scaled optimization space
    full_controls = maps_full.scaleVariables(full_controls_phys);
end

function schedule = updateScheduleWithControls(schedule, controls, maps)
    % Update schedule with new controls
    schedule = maps.updateSchedule(schedule, controls);
end

% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>