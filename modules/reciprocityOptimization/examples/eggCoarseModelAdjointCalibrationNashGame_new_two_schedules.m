%% Parameter tuning of a very coarse upscaling of the Egg model   
mrstModule add ad-core ad-blackoil deckformat agglom upscaling coarsegrid...
        mrst-gui ad-props incomp optimization test-suite linearsolvers 

%% Setup reference model
test    = TestCase('egg_wo', 'realization', 1);
problem = test.getPackedSimulationProblem();
problem.SimulatorSetup.model.useCNVConvergence = false;
problem.SimulatorSetup.model.nonlinearTolerance=1.0e-9;
problem.SimulatorSetup.model.toleranceMB =1.0e-9;
problem.SimulatorSetup.model.FacilityModel.toleranceWellBHP = 1000;
problem.SimulatorSetup.model.FacilityModel.toleranceWellRate = 5.0e-8;

% Focus on first 60 steps
schedule = problem.SimulatorSetup.schedule;
schedule.step.val     = schedule.step.val(1:60);
schedule.step.control = schedule.step.control(1:60);
Wref     = schedule.control.W;
dt       = schedule.step.val;
nstep    = numel(schedule.step.val);

% Create training schedule with varying controls for Player 1 (rate‑dominant)
perturbStep1 = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
rng(0)
schedule1 = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .4, 'perturbStep', perturbStep1);
problem.SimulatorSetup.schedule = schedule1;

% Create a DIFFERENT training schedule for Player 2 (pressure‑dominant)
% Use different random seed and different perturbation pattern
rng(42) % Different seed for different pattern
% Simulate reference for Player 1
problem.SimulatorSetup.schedule = schedule1;
simulatePackedProblem(problem);
[wsRef1, statesRef1] = getPackedSimulatorOutput(problem);
modelRef = problem.SimulatorSetup.model;

% Simulate reference for Player 2 (using schedule2 directly)
problem2 = test.getPackedSimulationProblem('Name', 'egg_wo_new2');
problem2.SimulatorSetup.model = modelRef; % Use same model settings
problem2.SimulatorSetup.model.useCNVConvergence = false;
problem2.SimulatorSetup.model.nonlinearTolerance=1.0e-9;
problem2.SimulatorSetup.model.toleranceMB =1.0e-9;
problem2.SimulatorSetup.model.FacilityModel.toleranceWellBHP = 1000;
problem2.SimulatorSetup.model.FacilityModel.toleranceWellRate = 5.0e-8;
problem2.SimulatorSetup.schedule = schedule;
simulatePackedProblem(problem2);
[wsRef2_2, statesRef2_2] = getPackedSimulatorOutput(problem2);
schedule2 = schedule;%flipWellControlsAndAddBC(schedule, statesRef2_2, []);
Wref2     = schedule2.control.W;
%perturbStep2 = [ones(1,8), round(.5+(2:(nstep-7))/2)]'; % Different pattern
M = 60;   % desired number of control steps (change as needed)
blockSize = ceil(nstep / M);
perturbStep2 = ceil((1:nstep) / blockSize)';
% Use different perturbation factors for pressure and rate
schedule2 = perturbedSimpleSchedule(dt, 'W', Wref2, ...
    'pressureFac', .1, 'rateFac', .3, 'perturbStep', perturbStep2);
problem2  = packSimulationProblem(problem.SimulatorSetup.state0, modelRef, schedule2, 'egg_wo_new2_2');
problem2.SimulatorSetup.model = modelRef; % Use same model settings
problem2.SimulatorSetup.model.useCNVConvergence = false;
problem2.SimulatorSetup.model.nonlinearTolerance=1.0e-9;
problem2.SimulatorSetup.model.toleranceMB =1.0e-9;
problem2.SimulatorSetup.model.FacilityModel.toleranceWellBHP = 1000;
problem2.SimulatorSetup.model.FacilityModel.toleranceWellRate = 5.0e-8;
simulatePackedProblem(problem2);
[wsRef2, statesRef2] = getPackedSimulatorOutput(problem2);
NashAlgorithm = true;

% Define playerChoice here at the top level for both Nash and non-Nash cases
if ~NashAlgorithm
    playerChoice = 2; % Set to 1 for Player 1 or 2 for Player 2
    fprintf('Non-Nash optimization selected. Optimizing %s\n', ...
        getPlayerName(playerChoice));
else
    playerChoice = 0; % Not used in Nash case, but defined to avoid errors
    fprintf('Nash optimization selected. Optimising both players with independent schedules.\n');
end

%% Create coarse models for both problems
% Common grid partitioning
blockIx = partitionUI(modelRef.G, [12, 12, 1]);
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);

% Coarse model for Player 1 (rate‑dominant)
modelCoarse1 = upscaleModelTPFA(modelRef, blockIx);
modelCoarse1.AutoDiffBackend = AutoDiffBackend();
pts1 = modelCoarse1.fluid.krPts;
scaling1 = {'SWL', pts1.w(1), 'SWCR', pts1.w(2), 'SWU', pts1.w(3), ...
           'SOWCR', pts1.ow(2), 'KRW', pts1.w(4), 'KRO', pts1.ow(4)};
modelCoarse1 = imposeRelpermScaling(modelCoarse1, scaling1{:});
modelCoarse1.toleranceCNV = 1e-6;

% Coarse model for Player 2 (pressure‑dominant) – using same upscaling method but with schedule2
modelCoarse2 = upscaleModelTPFA(modelRef, blockIx);
modelCoarse2.AutoDiffBackend = AutoDiffBackend();
pts2 = modelCoarse2.fluid.krPts;
scaling2 = {'SWL', pts2.w(1), 'SWCR', pts2.w(2), 'SWU', pts2.w(3), ...
           'SOWCR', pts2.ow(2), 'KRW', pts2.w(4), 'KRO', pts2.ow(4)};
modelCoarse2 = imposeRelpermScaling(modelCoarse2, scaling2{:});
modelCoarse2.toleranceCNV = 1e-6;

%% Set up initial states and schedules for coarse models
% Player 1 (rate‑dominant) problem
stateCoarse0_1 = upscaleState(modelCoarse1, modelRef, test.state0);
scheduleCoarse1 = upscaleSchedule(modelCoarse1, schedule1, 'wellUpscaleMethod', 'sum');
problem_coarse1 = packSimulationProblem(stateCoarse0_1, modelCoarse1, scheduleCoarse1, 'coarse_Player1');

% Player 2 (pressure‑dominant) problem – using schedule2 directly
stateCoarse0_2 = upscaleState(modelCoarse2, modelRef, test.state0);
scheduleCoarse2 = upscaleSchedule(modelCoarse2, schedule2, 'wellUpscaleMethod', 'sum');
problem_coarse2 = packSimulationProblem(stateCoarse0_2, modelCoarse2, scheduleCoarse2, 'coarse_Player2');

%% Simulate initial coarse models
[wsCoarse1, statesCoarse1] = simulateScheduleAD(stateCoarse0_1, modelCoarse1, scheduleCoarse1);
[wsCoarse2, statesCoarse2] = simulateScheduleAD(stateCoarse0_2, modelCoarse2, scheduleCoarse2);

% plot initial results
if NashAlgorithm
    plotWellSols({wsRef1, wsCoarse1, wsRef2, wsCoarse2}, ...
        {schedule1.step.val, scheduleCoarse1.step.val, schedule2.step.val, scheduleCoarse2.step.val},...
        'datasetnames',{'fine Player 1 (rate‑dominant)','coarse Player 1 (rate‑dominant)', ...
                        'fine Player 2 (pressure‑dominant)','coarse Player 2 (pressure‑dominant)'});
elseif playerChoice ==1
    plotWellSols({wsRef1, wsCoarse1}, ...
    {schedule1.step.val, scheduleCoarse1.step.val},...
    'datasetnames',{'fine scale model','initial upscaled Player 1 (rate‑dominant) model'});
elseif playerChoice == 2
    plotWellSols({wsRef2, wsCoarse2}, ...
    {schedule2.step.val, scheduleCoarse2.step.val},...
    'datasetnames',{'fine scale model','initial upscaled Player 2 (pressure‑dominant) model'});
end

%% Define parameter configurations for both problems
pv1 = modelCoarse1.operators.pv;
pv2 = modelCoarse2.operators.pv;

if NashAlgorithm
    % Parameters for Player 1 (rate‑dominant)
    config_player1 = {...
        %name           include    scaling              boxlims   relativeLimits
        'porevolume',       1,     'linear',    [.005*pv1, 2.5*pv1],              []
        'conntrans',        1,        'log',                  [],     [1e-3, 1e3]
        'transmissibility', 1,        'log',                  [],     [1e-3, 1e3]
        'swl',              1,     'linear',             [0, .5],              []
        'swcr',             1,     'linear',             [0, .5],              []
        'swu',              1,     'linear',             [.5, 1],              []
        'sowcr',            1,     'linear',             [0, .5],              []
        'krw',              1,     'linear',           [.25, 2.5],              []
        'kro',              1,     'linear',           [.25, 2.5],              []
        'sw',               1,     'linear',            [.005, 0.5],            []
        'pressure',         1,     'log',           [.895.*stateCoarse0_1.pressure, 1.5.*stateCoarse0_1.pressure], []
        };

    % Parameters for Player 2 (pressure‑dominant) – can use different parameters if desired
    config_player2 = {...
        %name           include    scaling              boxlims   relativeLimits
        'porevolume',       1,     'linear', [.005*pv2, 2.5*pv2],        []
        'conntrans',        1,        'log',                  [],     [1e-3, 1e3]
        'transmissibility', 1,        'log',                  [],     [1e-3, 1e3]
        'swl',              1,     'linear',             [0, .5],              []
        'swcr',             1,     'linear',             [0, .5],              []
        'swu',              1,     'linear',             [.5, 1],              []
        'sowcr',            1,     'linear',             [0, .5],              []
        'krw',              1,     'linear',           [.25, 2.5],              []
        'kro',              1,     'linear',           [.25, 2.5],              []
        'sw',               1,     'linear',           [.005, 0.5],             []
        'pressure',         1,     'log',           [.895.*stateCoarse0_2.pressure, 1.5.*stateCoarse0_2.pressure], []
        };
else
    % Non-Nash case – configure based on player choice
    if playerChoice == 1
        config_player1 = {...
            'porevolume',       1,     'linear',    [.005*pv1, 2.5*pv1],              []
            'conntrans',        1,        'log',                  [],     [1e-3, 1e3]
            'transmissibility', 1,        'log',                  [],     [1e-3, 1e3]
            'swl',              1,     'linear',             [0, .5],              []
            'swcr',             1,     'linear',             [0, .5],              []
            'swu',              1,     'linear',             [.5, 1],              []
            'sowcr',            1,     'linear',             [0, .5],              []
            'krw',              1,     'linear',           [.25, 2.5],              []
            'kro',              1,     'linear',           [.25, 2.5],              []
            'sw',               1,     'linear',            [.005, 0.5],            []
            'pressure',         1,     'log',           [.895.*stateCoarse0_1.pressure, 1.5.*stateCoarse0_1.pressure], []
            };
        config_player2 = []; % Not used
    else
        config_player2 = {...
            'porevolume',       1,     'linear',    [.005*pv2, 2.5*pv2],              []
            'conntrans',        1,        'log',                  [],     [1e-3, 1e3]
            'transmissibility', 1,        'log',                  [],     [1e-3, 1e3]
            'swl',              1,     'linear',             [0, .5],              []
            'swcr',             1,     'linear',             [0, .5],              []
            'swu',              1,     'linear',             [.5, 1],              []
            'sowcr',            1,     'linear',             [0, .5],              []
            'krw',              1,     'linear',           [.25, 2.5],              []
            'kro',              1,     'linear',           [.25, 2.5],              []
            'sw',               1,     'linear',            [.005, 0.5],            []
            'pressure',         1,     'log',           [.895.*stateCoarse0_2.pressure, 1.5.*stateCoarse0_2.pressure], []
            };
        config_player1 = []; % Not used
    end
end

%% Initialize parameters for both problems
parameters_player1 = [];
parameters_player2 = [];

if NashAlgorithm
    % Nash case: Initialize parameters for both players with their respective schedules
    for k = 1:size(config_player1,1)
        if config_player1{k, 2} == 0, continue, end
        parameters_player1 = addParameter(parameters_player1, struct('model', modelCoarse1, 'schedule', scheduleCoarse1, 'state0', stateCoarse0_1), ...
            'name', config_player1{k,1}, 'scaling', config_player1{k,3}, ...
            'boxLims', config_player1{k,4}, 'relativeLimits', config_player1{k,5});
    end
    
    for k = 1:size(config_player2,1)
        if config_player2{k, 2} == 0, continue, end
        parameters_player2 = addParameter(parameters_player2, struct('model', modelCoarse2, 'schedule', scheduleCoarse2, 'state0', stateCoarse0_2), ...
            'name', config_player2{k,1}, 'scaling', config_player2{k,3}, ...
            'boxLims', config_player2{k,4}, 'relativeLimits', config_player2{k,5});
    end
    
else
    % Non-Nash case: Only initialize parameters for the chosen player
    if playerChoice == 1
        for k = 1:size(config_player1,1)
            if config_player1{k, 2} == 0, continue, end
            parameters_player1 = addParameter(parameters_player1, struct('model', modelCoarse1, 'schedule', scheduleCoarse1, 'state0', stateCoarse0_1), ...
                'name', config_player1{k,1}, 'scaling', config_player1{k,3}, ...
                'boxLims', config_player1{k,4}, 'relativeLimits', config_player1{k,5});
        end
        parameters_player2 = [];
    else
        for k = 1:size(config_player2,1)
            if config_player2{k, 2} == 0, continue, end
            parameters_player2 = addParameter(parameters_player2, struct('model', modelCoarse2, 'schedule', scheduleCoarse2, 'state0', stateCoarse0_2), ...
                'name', config_player2{k,1}, 'scaling', config_player2{k,3}, ...
                'boxLims', config_player2{k,4}, 'relativeLimits', config_player2{k,5});
        end
        parameters_player1 = [];
    end
end

%% Define mismatch functions for each twin
weighting1  = {'WaterRateWeight', 1.0/(150/day), ...
             'OilRateWeight',   1.0/(80/day), ...
             'BHPWeight',      1/(20*barsa)};
weighting2 = weighting1;
% Compute weights from the respective reference well solutions
[waterWeight1, oilWeight1, bhpWeight1] = computeWeightsFromWellSol(wsRef1);
[waterWeight2, oilWeight2, bhpWeight2] = computeWeightsFromWellSol(wsRef2);

% Build weighting cell arrays for each problem
weighting1 = {'WaterRateWeight', waterWeight1, ...
              'OilRateWeight',   oilWeight1,   ...
              'BHPWeight',       bhpWeight1};

weighting2 = {'WaterRateWeight', waterWeight2, ...
              'OilRateWeight',   oilWeight2,   ...
              'BHPWeight',       bhpWeight2};
% Create objective functions that match unitBoxLM signature
player1_solver ='unitboxlm';
player2_solver = 'unitboxlm';

if NashAlgorithm || playerChoice ==1
    if strcmp(player1_solver, 'unitboxlm')
            % Objective for Player 1 (rate‑dominant): Match reference data from schedule1
            mismatchFn1 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
                matchObservedOW(model, states, schedule, states_ref,...
                'computePartials', compDer, 'tstep', tstep, weighting1{:},...
                'state', state, 'from_states', false, 'mismatchSum', false);
            obj_fn1 = @(p) evaluateMatchSummands(p, mismatchFn1, ...
                problem_coarse1.SimulatorSetup, parameters_player1, statesRef1);
    else
            mismatchFn1 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
                matchObservedOW(model, states, schedule, states_ref,...
                'computePartials', compDer, 'tstep', tstep, weighting1{:},...
                'state', state, 'from_states', false);
            obj_fn1 = @(p) evaluateMatch(p, mismatchFn1, ...
                problem_coarse1.SimulatorSetup, parameters_player1, statesRef1);
    end
end

if NashAlgorithm || playerChoice ==2
    if strcmp(player2_solver, 'unitboxlm')
        % Objective for Player 2 (pressure‑dominant): Match reference data from schedule2
        mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
            matchObservedOW(model, states, schedule, states_ref,...
            'computePartials', compDer, 'tstep', tstep, weighting2{:},...
            'state', state, 'from_states', false, 'mismatchSum', false);
        obj_fn2 = @(p) evaluateMatchSummands(p, mismatchFn2, ...
            problem_coarse2.SimulatorSetup, parameters_player2, statesRef2);
    else
        mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
            matchObservedOW(model, states, schedule, states_ref,...
            'computePartials', compDer, 'tstep', tstep, weighting2{:},...
            'state', state, 'from_states', false);
        obj_fn2 = @(p) evaluateMatch(p, mismatchFn2, ...
            problem_coarse2.SimulatorSetup, parameters_player2, statesRef2);
    end
end

%% Create and run Nash optimization using NashOptimizationTwinProblem class

% Get initial parameter vectors
if NashAlgorithm || playerChoice ==1
    p0_player1 = getScaledParameterVector(problem_coarse1.SimulatorSetup, parameters_player1);
end

if NashAlgorithm || playerChoice ==2
    p0_player2 = getScaledParameterVector(problem_coarse2.SimulatorSetup, parameters_player2);
end

if NashAlgorithm
    % Create Nash optimization problem with independent schedules
    nashProblem = NashOptimizationTwinProblem(...
        'twinParameters', {parameters_player1, parameters_player2}, ...
        'twinProblems', {problem_coarse1, problem_coarse2}, ...
        'twinReferenceStates', {statesRef1, statesRef2}, ...
        'twinObjectives', {mismatchFn1, mismatchFn2}, ...
        'nashParameters', struct('maxIterations', 15, ...
        'convergenceTol', 1e-4, ...
        'mixingParameter', 0.65, ...
        'regularizationParam', [1.5e-5,1.5e-4], ...
        'maxInnerIterations', 10, ...
        'verbose', true));

    % Run optimization
    initialGuess = {p0_player1, p0_player2};
    [u_opt, history, nashResults] = nashProblem.optimizeNashEquilibrium(initialGuess, ...
        'optimizer_1', 'unitboxlm', 'optimizer_2', 'unitboxlm', 'objectiveFunctions', {obj_fn1, obj_fn2});

    % Get optimal setups
    [setup_opt1, setup_opt2] = nashProblem.getOptimalSetups();

    % Plot convergence
    nashProblem.plotTwinConvergence();

    %% Simulate final tuned models
    [wellSols_opt1, states_opt1] = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, setup_opt1.schedule);
    [wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);

    % Plot results
    plotWellSols({wsRef1, wellSols_opt1}, ...
        repmat({schedule1.step.val}, 1, 2), ...
        'datasetnames', {'reference','coarse tuned Player 1 (rate‑dominant)'});

    plotWellSols({wsRef2, wellSols_opt2}, ...
        repmat({schedule2.step.val}, 1, 2), ...
        'datasetnames', {'reference','coarse tuned Player 2 (pressure‑dominant)'});

    plotWellSols({wsRef1, wsRef2, wellSols_opt1, wellSols_opt2}, ...
        {schedule1.step.val, schedule2.step.val, schedule1.step.val, schedule2.step.val}, ...
        'datasetnames', {'ref Player 1 (rate‑dominant)','ref Player 2 (pressure‑dominant)', ...
                         'tuned Player 1 (rate‑dominant)','tuned Player 2 (pressure‑dominant)'});
    
    %% Display final results
    fprintf('\n=== Optimization Results ===\n');
    fprintf('Final objective values:\n');
    fprintf('  Player 1 (rate‑dominant) problem: %.4e\n', history(end).twin{1}.val(end));
    fprintf('  Player 2 (pressure‑dominant) problem: %.4e\n', history(end).twin{2}.val(end));

    fprintf('\nParameter changes from initial to final:\n');
    param_change1 = norm(u_opt{1} - p0_player1) / norm(p0_player1);
    param_change2 = norm(u_opt{2} - p0_player2) / norm(p0_player2);
    fprintf('  Player 1 (rate‑dominant) parameters: %.2f%% change\n', param_change1 * 100);
    fprintf('  Player 2 (pressure‑dominant) parameters: %.2f%% change\n', param_change2 * 100);
    
    fprintf('\nSchedule information:\n');
    fprintf('  Player 1 (rate‑dominant) uses schedule1 with %d steps\n', numel(schedule1.step.val));
    fprintf('  Player 2 (pressure‑dominant) uses schedule2 with %d steps\n', numel(schedule2.step.val));
else
    % Non-Nash optimization for the chosen player
    if playerChoice == 1
        fprintf('Optimising Player 1 (rate‑dominant) problem with schedule1...\n');
        if strcmp(player1_solver,'unitboxlm')
            [v1, p_opt1, history1] = unitBoxLM(p0_player1, obj_fn1, 'maxIt', 20);
        else
            [v1, p_opt1, history1] = unitBoxBFGS(p0_player1, obj_fn1, 'objChangeTol', 1e-5, ...
                                      'maxIt', 25, 'logPlot', true);
        end
        setup_opt1 = updateSetupFromScaledParameters(problem_coarse1.SimulatorSetup, parameters_player1, p_opt1);
        [wellSols_opt1, states_opt1] = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, setup_opt1.schedule);
        % Plot results
        plotWellSols({wsRef1, wellSols_opt1}, ...
            repmat({schedule1.step.val}, 1, 2), ...
            'datasetnames', {'reference','coarse tuned Player 1 (rate‑dominant)'});

        %% Display final results
        fprintf('\n=== Optimization Results ===\n');
        fprintf('Final objective values:\n');
        fprintf('  Player 1 (rate‑dominant) problem: %.4e\n', history1(end).val(end));
        fprintf('\nParameter changes from initial to final:\n');
        param_change1 = norm(p_opt1 - p0_player1) / norm(p0_player1);
        fprintf('  Player 1 (rate‑dominant) parameters: %.2f%% change\n', param_change1 * 100);
    else
        fprintf('Optimising Player 2 (pressure‑dominant) problem with schedule2...\n');
        if strcmp(player2_solver,'unitboxlm')
            [v2, p_opt2, history2] = unitBoxLM(p0_player2, obj_fn2, 'maxIt', 20);
        else
            [v2, p_opt2, history2] = unitBoxBFGS(p0_player2, obj_fn2, 'objChangeTol', 1e-5, ...
                                      'maxIt', 25, 'logPlot', true);
        end
        setup_opt2 = updateSetupFromScaledParameters(problem_coarse2.SimulatorSetup, parameters_player2, p_opt2);
        [wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);
        % Plot results
        plotWellSols({wsRef2, wellSols_opt2}, ...
            repmat({schedule2.step.val}, 1, 2), ...
            'datasetnames', {'reference','coarse tuned Player 2 (pressure‑dominant)'});

        %% Display final results
        fprintf('\n=== Optimization Results ===\n');
        fprintf('Final objective values:\n');
        fprintf('  Player 2 (pressure‑dominant) problem: %.4e\n', history2(end).val(end));
        fprintf('\nParameter changes from initial to final:\n');
        param_change2 = norm(p_opt2 - p0_player2) / norm(p0_player2);
        fprintf('  Player 2 (pressure‑dominant) parameters: %.2f%% change\n', param_change2 * 100);
    end
end

%% Compare reference, initial coarse and optimized coarse for validation schedules
rng(100);
W = test.schedule.control.W;

% Create validation schedule for Player 1 (rate‑dominant) – different from training
s_valid1 = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .2, 'perturbStep', ones(numel(dt),1));

% Create validation schedule for Player 2 (pressure‑dominant) – different from training
s_valid2 = perturbedSimpleSchedule(dt, 'W', Wref2, ...
    'pressureFac', .1, 'rateFac', .2, 'perturbStep', ones(numel(dt),1));

% Simulate reference on validation schedules
[ws_valid1, states_valid1] = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_valid1, ...
    'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);
[ws_valid2, states_valid2] = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_valid2, ...
    'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);

if NashAlgorithm
    %% Nash Algorithm Case – Validate both players on their respective validation schedules
    
    % Player 1 (rate‑dominant) validation
    w_opt1 = setup_opt1.schedule.control(1).W;
    w_old1 = scheduleCoarse1.control(1).W;
    s_opt1_valid = simpleSchedule(dt, 'W', w_opt1);
    s_old1_valid = simpleSchedule(dt, 'W', w_old1);
    
    % Apply Player 1 validation controls
    for kw = 1:numel(Wref)
        s_opt1_valid.control.W(kw).val = s_valid1.control.W(kw).val;
        s_old1_valid.control.W(kw).val = s_valid1.control.W(kw).val;
    end
    
    % Player 2 (pressure‑dominant) validation
    w_opt2 = setup_opt2.schedule.control(1).W;
    w_old2 = scheduleCoarse2.control(1).W;
    s_opt2_valid = simpleSchedule(dt, 'W', w_opt2);
    s_old2_valid = simpleSchedule(dt, 'W', w_old2);
    
    % Apply Player 2 validation controls
    for kw = 1:numel(Wref)
        s_opt2_valid.control.W(kw).val = s_valid2.control.W(kw).val;
        s_old2_valid.control.W(kw).val = s_valid2.control.W(kw).val;
    end
    
    % Simulate predictions for both players on validation schedules
    ws_coarse_valid1_opt = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, s_opt1_valid);
    ws_coarse_valid1_old = simulateScheduleAD(stateCoarse0_1, modelCoarse1, s_old1_valid);
    
    ws_coarse_valid2_opt = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, s_opt2_valid);
    ws_coarse_valid2_old = simulateScheduleAD(stateCoarse0_2, modelCoarse2, s_old2_valid);
    
    % Plot results for both players on validation schedules
    figure('Position', [100, 100, 1200, 500]);
    subplot(1,2,1);
    plotWellSols({ws_valid1, ws_coarse_valid1_old, ws_coarse_valid1_opt}, ...
        repmat({schedule1.step.val}, 1, 3), ...
        'datasetnames', {'reference','coarse initial','coarse tuned'});
    title('Player 1 (rate‑dominant) – Validation with independent schedule');
    
    subplot(1,2,2);
    plotWellSols({ws_valid2, ws_coarse_valid2_old, ws_coarse_valid2_opt}, ...
        repmat({schedule2.step.val}, 1, 3), ...
        'datasetnames', {'reference','coarse initial','coarse tuned'});
    title('Player 2 (pressure‑dominant) – Validation with independent schedule');
    
else
    %% Non-Nash Case – Validate only the chosen player
    
    if playerChoice == 1
        % Player 1 (rate‑dominant) validation
        w_opt1 = setup_opt1.schedule.control(1).W;
        w_old1 = scheduleCoarse1.control(1).W;
        s_opt1_valid = simpleSchedule(dt, 'W', w_opt1);
        s_old1_valid = simpleSchedule(dt, 'W', w_old1);
        
        for kw = 1:numel(Wref)
            s_opt1_valid.control.W(kw).val = s_valid1.control.W(kw).val;
            s_old1_valid.control.W(kw).val = s_valid1.control.W(kw).val;
        end
        
        ws_coarse_valid1_opt = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, s_opt1_valid);
        ws_coarse_valid1_old = simulateScheduleAD(stateCoarse0_1, modelCoarse1, s_old1_valid);
        
        % Plot results for Player 1
        plotWellSols({ws_valid1, ws_coarse_valid1_old, ws_coarse_valid1_opt}, ...
            repmat({schedule1.step.val}, 1, 3), ...
            'datasetnames', {'reference','coarse initial','coarse tuned'});
        title('Player 1 (rate‑dominant) – Validation with independent schedule');
        
    else
        % Player 2 (pressure‑dominant) validation
        w_opt2 = setup_opt2.schedule.control(1).W;
        w_old2 = scheduleCoarse2.control(1).W;
        s_opt2_valid = simpleSchedule(dt, 'W', w_opt2);
        s_old2_valid = simpleSchedule(dt, 'W', w_old2);
        
        for kw = 1:numel(Wref)
            s_opt2_valid.control.W(kw).val = s_valid2.control.W(kw).val;
            s_old2_valid.control.W(kw).val = s_valid2.control.W(kw).val;
        end
        
        ws_coarse_valid2_opt = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, s_opt2_valid);
        ws_coarse_valid2_old = simulateScheduleAD(stateCoarse0_2, modelCoarse2, s_old2_valid);
        
        % Plot results for Player 2
        plotWellSols({ws_valid2, ws_coarse_valid2_old, ws_coarse_valid2_opt}, ...
            repmat({schedule2.step.val}, 1, 3), ...
            'datasetnames', {'reference','coarse initial','coarse tuned'});
        title('Player 2 (pressure‑dominant) – Validation with independent schedule');
    end
end

%% Store final models for later use
if NashAlgorithm
    final_models.nash = true;
    final_models.player1.model = setup_opt1.model;
    final_models.player1.state0 = setup_opt1.state0;
    final_models.player1.schedule = setup_opt1.schedule;
    final_models.player1.wellSols = wellSols_opt1;
    final_models.player1.states = states_opt1;
    final_models.player1.schedule_training = schedule1;
    final_models.player1.states_ref = statesRef1;
    
    final_models.player2.model = setup_opt2.model;
    final_models.player2.state0 = setup_opt2.state0;
    final_models.player2.schedule = setup_opt2.schedule;
    final_models.player2.wellSols = wellSols_opt2;
    final_models.player2.states = states_opt2;
    final_models.player2.schedule_training = schedule2;
    final_models.player2.states_ref = statesRef2;
    
    fprintf('\n=== Final Models Stored ===\n');
    fprintf('Both Player 1 (rate‑dominant) and Player 2 (pressure‑dominant) models are available\n');
    fprintf('  Player 1 trained on schedule1, validated on independent schedule\n');
    fprintf('  Player 2 trained on schedule2, validated on independent schedule\n');
    
else
    final_models.nash = false;
    if playerChoice == 1
        final_models.player1.model = setup_opt1.model;
        final_models.player1.state0 = setup_opt1.state0;
        final_models.player1.schedule = setup_opt1.schedule;
        final_models.player1.wellSols = wellSols_opt1;
        final_models.player1.states = states_opt1;
        final_models.player1.schedule_training = schedule1;
        final_models.player1.states_ref = statesRef1;
        fprintf('\n=== Final Model Stored ===\n');
        fprintf('Player 1 (rate‑dominant) model trained on schedule1\n');
    else
        final_models.player2.model = setup_opt2.model;
        final_models.player2.state0 = setup_opt2.state0;
        final_models.player2.schedule = setup_opt2.schedule;
        final_models.player2.wellSols = wellSols_opt2;
        final_models.player2.states = states_opt2;
        final_models.player2.schedule_training = schedule2;
        final_models.player2.states_ref = statesRef2;
        fprintf('\n=== Final Model Stored ===\n');
        fprintf('Player 2 (pressure‑dominant) model trained on schedule2\n');
    end
end

function [waterWeight, oilWeight, bhpWeight] = computeWeightsFromWellSol(wsRef)
    % wsRef: cell array of well solutions (one per timestep)
    % Returns weights = 1 / (max - min) for each quantity
    
    qWs_all = [];
    qOs_all = [];
    bhp_all = [];
    
    for step = 1:numel(wsRef)
        ws_step = wsRef{step};
        for w = 1:numel(ws_step)
            qWs_all = [qWs_all; ws_step(w).qWs];
            qOs_all = [qOs_all; ws_step(w).qOs];
            bhp_all = [bhp_all; ws_step(w).bhp];
        end
    end
    
    % Water rate weight
    rW = max(qWs_all) - min(qWs_all);
    if rW == 0, rW = 1; end
    waterWeight = 1 / rW;
    
    % Oil rate weight
    rO = max(qOs_all) - min(qOs_all);
    if rO == 0, rO = 1; end
    oilWeight = 1 / rO;
    
    % BHP weight
    rB = max(bhp_all) - min(bhp_all);
    if rB == 0, rB = 1; end
    bhpWeight = 1 / rB;
end

function name = getPlayerName(playerChoice)
    if playerChoice == 1
        name = 'Player 1 (rate‑dominant perturbation)';
    else
        name = 'Player 2 (pressure‑dominant perturbation)';
    end
end