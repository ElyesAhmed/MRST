mrstModule add ad-core ad-blackoil ad-props optimization spe10
setupModel2D

% Create model
model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
state0 = initResSol(G, 200*barsa, [0, 1]);

problem = packSimulationProblem(state0, model, schedule, 'tmp', 'Name', 'model2D', 'Description', 'model2d');

numRealizations = 1;
problems = cell(1, numRealizations);
for k = 1:numRealizations
    problems{k} = problem;
    problems{k}.seed = k;
end
samples = struct('problem', {problems}, 'num', numRealizations);

% Define objective functions
ro = 60;      % oil revenue/price ($/stb)
rwp = 6;      % water production handling costs ($/stb)
rwi = 6;      % water injection cost ($/stb) 

% Discounted NPV objective
objectiveOptsDiscounted = {'OilPrice', ro, ...
                          'WaterInjectionCost', rwi, ...
                          'WaterProductionCost', rwp, ...
                          'DiscountFactor', 0.25};

% Undiscounted NPV objective                          
objectiveOptsUndiscounted = {'OilPrice', ro, ...
                            'WaterInjectionCost', rwi, ...
                            'WaterProductionCost', rwp, ...
                            'DiscountFactor', 0};

objFunDiscounted = @(model, states, schedule, varargin) ...
    NPVOW(model, states, schedule, varargin{:}, objectiveOptsDiscounted{:});

objFunUndiscounted = @(model, states, schedule, varargin) ...
    NPVOW(model, states, schedule, varargin{:}, objectiveOptsUndiscounted{:});

% Objective scaling
objectiveScaling = 3.2e7;

% Create objective structures
objective1 = struct('function', objFunDiscounted, 'scaling', objectiveScaling);
objective2 = struct('function', objFunUndiscounted, 'scaling', objectiveScaling);

% Your bounds setup
bounds = [];
for j = 1:length(schedule.control)    
    W = schedule.control(j).W;    
    for i = 1:length(W)
        if W(i).sign > 0  % Injector
            ref_rate = W(i).val;            
            step_bounds = [0, ref_rate*3];
            bounds = [bounds; step_bounds];
        end
    end
    for i = 1:length(W)
        if W(i).sign < 0  % Producer
            step_bounds = [100, 250] * barsa;
            bounds = [bounds; step_bounds];
        end
    end
end

% Setup mappings
bnds = processBounds(W, ...
    'rate(inj)', bounds(1,:), ...
    'rate(inj)', bounds(2,:), ...
    'bhp(prod)', bounds(3,:), ...
    'bhp(prod)', bounds(4,:),  ...
    'rate(inj)', bounds(5,:), ...
    'rate(inj)', bounds(6,:), ...
    'bhp(prod)', bounds(7,:), ...
    'bhp(prod)', bounds(8,:), ...
    'rate(inj)', bounds(9,:), ...
    'rate(inj)', bounds(10,:), ...
    'bhp(prod)', bounds(11,:), ...
    'bhp(prod)', bounds(12,:), ...
    'rate(inj)', bounds(13,:), ...
    'rate(inj)', bounds(14,:), ...
    'bhp(prod)', bounds(15,:), ...
    'bhp(prod)', bounds(16,:));

[fullMaps, fullUInit] = setupSimulationControlMappings(schedule, bnds);
[fullMaps, fullUInit] = reduceControlMappings(fullMaps, fullUInit, 'exludeLims', true);

% Define control types for each objective
controlTypes = {{'bhp', 'lrat'}, {'rate', 'rate'}};

% Create Nash optimization problem
nashProblem = NashOptimizationProblem(samples, ...
    'name', 'Nash_NPV_Optimization', ...
    'objectives', {objective1, objective2}, ...
    'maps', fullMaps, ...          % Full mapping
    'u_init', fullUInit, ...       % Full initial controls
    'bounds', bounds, ...          % Full bounds
    'controlTypes', controlTypes, ...  % Control types for each objective
    'setupType', 'simulation', ...
    'nashParameters', struct('maxIterations', 100, ...
                            'convergenceTol', 1e-4, ...
                            'mixingParameter', 0.6, ...
                            'verbose', true), ...
    'background', true);

% Run Nash equilibrium optimization - pass the full initial controls
[u_opt, history, nashResults] = nashProblem.optimizeNashEquilibrium(problems{1}, ...
    'maxIterations', 100, ...
    'convergenceTol', 1e-4, ...
    'mixingParameter', 0.6, ...
     'stepInit', 0.1, 'maxIt', 20, 'regularizationParam', 0.025);

% Display results
fprintf('\n=== Nash Equilibrium Results ===\n');
fprintf('Final controls: %s\n', mat2str(u_opt, 3));
fprintf('Final objective values:\n');
fprintf('  Discounted NPV: $%.2f\n', nashResults.objectiveValues(end, 1));
fprintf('  Undiscounted NPV: $%.2f\n', nashResults.objectiveValues(end, 2));

% Plot convergence using the new method
nashProblem.plotNashConvergence(nashResults);

% Simulate final optimal schedule
finalProblem = nashProblem.getSampleProblem(1);
finalProblem = nashProblem.updateProblemFun(finalProblem, u_opt);

[wellSols, states] = simulateScheduleAD(state0, model, finalProblem.SimulatorSetup.schedule);

% Calculate final NPV values
final_npv_disc = objFunDiscounted(model, states, finalProblem.SimulatorSetup.schedule);
final_npv_undisc = objFunUndiscounted(model, states, finalProblem.SimulatorSetup.schedule);

fprintf('\nFinal simulation results:\n');
fprintf('  Discounted NPV: $%.2f\n', sum([final_npv_disc{:}]));
fprintf('  Undiscounted NPV: $%.2f\n', sum([final_npv_undisc{:}]));

% Get objective values for specific iteration
iteration = 5;
obj_values = nashResults.objectiveValues(iteration, :);

% Get controls for specific iteration  
controls = nashResults.controls{iteration};

% Plot individual objective progress (inherited from parent)
nashProblem.plotObjectiveValues();

% Reset and restart optimization if needed
nashProblem.reset('prompt', false);
[u_opt_restart, history_restart] = nashProblem.optimizeNashEquilibrium(u_opt, ...
    'maxIterations', 5, 'convergenceTol', 1e-5);