%% optimizeModel1D - NPV Optimization for 1D Reservoir Model
% Optimize well controls for maximum NPV using different discount factors
%close all;
clear all;
mrstModule add ad-core ad-blackoil ad-props mrst-gui optimization ensemble

%% Setup optimization problem
% Setup base problems
numRealizations = 1;
problems = cell(1, numRealizations);
for k = 1:numRealizations
    problems{k} = setupFnSimpleProblem1D(k);
    problems{k}.seed = k;
end
samples = struct('problem', {problems}, 'num', numRealizations);

% Get schedule and model from first problem
schedule = problems{1}.SimulatorSetup.schedule;
model = problems{1}.SimulatorSetup.model;
state0 = problems{1}.SimulatorSetup.state0;
G = model.G;

% Setup target and limit bounds
W = schedule.control(1).W;
irate = W(1).val;
bnds = processBounds(W, ...
    'rate(inj)', [.1*irate, 3*irate], ...     % allowed target rates
    'bhp(inj)', [120, 500]*barsa, ...         % allowed upper bhp-limits
    'bhp(prod)', [50, 120]*barsa, ...         % allowed target bhp
    'lrat(prod)', [-3*irate, -.1*irate]);     % allowed lower lrat-limit

[maps, u_init] = setupSimulationControlMappings(schedule, bnds);

% Objective scaling
objectiveScaling = .9 * 40 * sum(G.cells.volumes) / stb;

%% Define objective functions
objectiveOptsDiscounted = {'OilPrice', 40/stb, ...
                          'WaterInjectionCost', 5/stb, ...
                          'WaterProductionCost', 4/stb, ...
                          'DiscountFactor', 0.5};

objectiveOptsUndiscounted = {'OilPrice', 40/stb, ...
                            'WaterInjectionCost', 5/stb, ...
                            'WaterProductionCost', 4/stb, ...
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
[maps_disc, u_disc_init] = reduceControlMappings(maps, u_init, 'type', {'rate', 'lrat'});
[maps_undisc, u_undisc_init] = reduceControlMappings(maps, u_init, 'type', {'bhp', 'bhp'});

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
u_disc_current = u_disc_init/u_disc_init;
u_undisc_current = u_undisc_init/u_undisc_init;

% Nash optimization parameters
maxNashIterations = 10;
convergence_tolerance = 1e-4;
nashResults = struct();

% Initialize problem for Nash iterations
problem_current = problems{1};

for nashIter = 1:maxNashIterations
    fprintf('\n--- Nash Iteration %d/%d ---\n', nashIter, maxNashIterations);
    
    % Phase 1: Optimize injection controls for discounted NPV
    fprintf('Optimizing injection controls (discounted NPV)...\n');
    p_discounted.reset('prompt', false)
    [u_disc_opt, history_disc] = p_discounted.optimize(...
        problem_current, 'stepInit', 0.1, 'maxIt', 20);
    
    % Update problem with optimized injection controls

    u_disc_current =  p_discounted.unscaleVariables(history_disc.u{1});
    problem_current = p_discounted.updateProblemFun(problem_current, 0.25*u_disc_opt+0.75*u_disc_current);    
    % Phase 2: Optimize production controls for undiscounted NPV
    fprintf('Optimizing production controls (undiscounted NPV)...\n');
    p_undiscounted.reset('prompt', false)
    [u_undisc_opt, history_undisc] = p_undiscounted.optimize(...
        problem_current, 'stepInit', 0.1, 'maxIt', 20);
    
    % Update problem with optimized production controls

    u_undisc_current =  p_undiscounted.unscaleVariables(history_undisc.u{1});

    problem_current = p_undiscounted.updateProblemFun(problem_current, 0.25*u_undisc_opt+0.75*u_undisc_current);    
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