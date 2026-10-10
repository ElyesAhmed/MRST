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

% Create training schedule with varying controls
perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
rng(0)
schedule = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .4, 'perturbStep', perturbStep);
problem.SimulatorSetup.schedule = schedule;

% Simulate
simulatePackedProblem(problem);
[wsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef = problem.SimulatorSetup.model;
NashAlgorithm = true;

% Define playerChoice here at the top level for both Nash and non-Nash cases
if ~NashAlgorithm
    playerChoice = 2; % Set to 1 for RATE-BHP (Player 1) or 2 for BHP-RATE (Player 2)
    fprintf('Non-Nash optimization selected. Optimizing Player %d\n', playerChoice);
else
    playerChoice = 0; % Not used in Nash case, but defined to avoid errors
    fprintf('Nash optimization selected. Optimizing both players.\n');
end

if NashAlgorithm || playerChoice ==2
    % Create schedule with flipped well controls
    schedule_new = flipWellControlsAndAddBC(schedule, statesRef, []);
    problem_new  = packSimulationProblem(problem.SimulatorSetup.state0, modelRef, schedule_new, 'egg_wo_new');
    problem_new.SimulatorSetup.model = modelRef; % Use same model settings
    simulatePackedProblem(problem_new);
    [wsRef_new, statesRef_new] = getPackedSimulatorOutput(problem_new);
end
%% Create coarse models for both problems
% Common grid partitioning
blockIx = partitionUI(modelRef.G, [12, 12, 1]);
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);

% Coarse model for RATE-Inj-BHP-Prod problem
modelCoarse = upscaleModelTPFA(modelRef, blockIx);
modelCoarse.AutoDiffBackend = AutoDiffBackend();
pts = modelCoarse.fluid.krPts;
scaling = {'SWL', pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
           'SOWCR', pts.ow(2), 'KRW', pts.w(4), 'KRO', pts.ow(4)};
modelCoarse = imposeRelpermScaling(modelCoarse, scaling{:});
modelCoarse.toleranceCNV = 1e-6;

if NashAlgorithm || playerChoice ==2
    % Coarse model for BHP-Inj-RATE-Prod problem
    modelCoarse_new = upscaleModelTPFA(modelRef, blockIx);
    modelCoarse_new.AutoDiffBackend = AutoDiffBackend();
    pts_new = modelCoarse_new.fluid.krPts;
    scaling_new = {'SWL', pts_new.w(1), 'SWCR', pts_new.w(2), 'SWU', pts_new.w(3), ...
        'SOWCR', pts_new.ow(2), 'KRW', pts_new.w(4), 'KRO', pts_new.ow(4)};
    modelCoarse_new = imposeRelpermScaling(modelCoarse_new, scaling_new{:});
    modelCoarse_new.toleranceCNV = 1e-6;
end
%% Set up initial states and schedules for coarse models
% RATE-Inj-BHP-Prod problem
stateCoarse0 = upscaleState(modelCoarse, modelRef, test.state0);
scheduleCoarse = upscaleSchedule(modelCoarse, schedule, 'wellUpscaleMethod', 'sum');
problem_coarse = packSimulationProblem(stateCoarse0, modelCoarse, scheduleCoarse, 'coarse_RATE_BHP');

if NashAlgorithm || playerChoice ==2
    % BHP-Inj-RATE-Prod problem
    stateCoarse0_new = upscaleState(modelCoarse_new, modelRef, test.state0);
    scheduleCoarse_new = upscaleSchedule(modelCoarse_new, schedule_new, 'wellUpscaleMethod', 'sum');
    problem_coarse_new = packSimulationProblem(stateCoarse0_new, modelCoarse_new, scheduleCoarse_new, 'coarse_BHP_RATE');
end
%% Simulate initial coarse models
[wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);
if NashAlgorithm || playerChoice ==2
    [wsCoarse_new, statesCoarse_new] = simulateScheduleAD(stateCoarse0_new, modelCoarse_new, scheduleCoarse_new);
end
% plot initial results
if NashAlgorithm
plotWellSols({wsRef, wsCoarse, wsCoarse_new}, ...
    {schedule.step.val, scheduleCoarse.step.val, scheduleCoarse_new.step.val},...
    'datasetnames',{'fine scale model','initial upscaled RATE-BHP model','initial upscaled BHP-RATE model'});
elseif playerChoice ==1
    plotWellSols({wsRef, wsCoarse}, ...
    {schedule.step.val, scheduleCoarse.step.val},...
    'datasetnames',{'fine scale model','initial upscaled RATE-BHP model'});
elseif playerChoice == 2
    plotWellSols({wsRef, wsCoarse_new}, ...
    {schedule.step.val, scheduleCoarse_new.step.val},...
    'datasetnames',{'fine scale model','initial upscaled BHP-RATE model'});
end
%% Define parameter configurations for both problems
pv = modelCoarse.operators.pv;
if NashAlgorithm || playerChoice ==2
    pv_new = modelCoarse_new.operators.pv;
end

if NashAlgorithm
    % Parameters for RATE-BHP problem (focus on transmissibility)
    config_player1 = {...
        %name           include    scaling              boxlims   relativeLimits
        'porevolume',       0,     'linear',    [.005*pv, 2.5*pv],              []
        'conntrans',        0,        'log',                  [],     [1e-3, 1e3]
        'transmissibility', 1,        'log',                  [],     [1e-3, 1e3]
        'swl',              1,     'linear',             [0, .5],              []
        'swcr',             1,     'linear',             [0, .5],              []
        'swu',              0,     'linear',             [.5, 1],              []
        'sowcr',            0,     'linear',             [0, .5],              []
        'krw',              0,     'linear',           [.25, 2.5],              []
        'kro',              1,     'linear',           [.25, 2.5],              []
        'sw',               1,     'linear',            [.005, 0.5],            []
        'pressure',         0,     'log',           [.895.*stateCoarse0.pressure, 1.5.*stateCoarse0.pressure], []
        };

    % Parameters for BHP-RATE problem (focus on relperm and pore volume)
    config_player2 = {...
        %name           include    scaling              boxlims   relativeLimits
        'porevolume',       1,     'linear', [.005*pv_new, 2.5*pv_new],        []
        'conntrans',        1,        'log',                  [],     [1e-3, 1e3]
        'transmissibility', 0,        'log',                  [],     [1e-3, 1e3]
        'swl',              0,     'linear',             [0, .5],              []
        'swcr',             0,     'linear',             [0, .5],              []
        'swu',              1,     'linear',             [.5, 1],              []
        'sowcr',            1,     'linear',             [0, .5],              []
        'krw',              1,     'linear',           [.25, 2.5],              []
        'kro',              0,     'linear',           [.25, 2.5],              []
        'sw',               0,     'linear',           [.005, 0.5],             []
        'pressure',         1,     'log',           [.895.*stateCoarse0.pressure, 1.5.*stateCoarse0.pressure], []
        };
else

    % Player 1 (RATE-BHP) - comprehensive parameters
    config_player1 = {...
        %name           include    scaling              boxlims   relativeLimits
        'porevolume',       1,     'linear',    [.005*pv, 2.5*pv],              []
        'conntrans',        0,        'log',                  [],     [1e-3, 1e3]
        'transmissibility', 0,        'log',                  [],     [1e-3, 1e3]
        'swl',              0,     'linear',             [0, .5],              []
        'swcr',             0,     'linear',             [0, .5],              []
        'swu',              1,     'linear',             [.5, 1],              []
        'sowcr',            1,     'linear',             [0, .5],              []
        'krw',              0,     'linear',           [.25, 2.5],              []
        'kro',              1,     'linear',           [.25, 2.5],              []
        'sw',               0,     'linear',            [.005, 0.5],            []
        'pressure',         1,     'log',           [.895.*stateCoarse0.pressure, 1.5.*stateCoarse0.pressure], []
        };
    config_player2 = config_player1;

end

%% Initialize parameters for both problems
parameters_player1 = [];
parameters_player2 = [];

if NashAlgorithm
    % Nash case: Initialize parameters for both players
    for k = 1:size(config_player1,1)
        if config_player1{k, 2} == 0, continue, end
        parameters_player1 = addParameter(parameters_player1, struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0), ...
            'name', config_player1{k,1}, 'scaling', config_player1{k,3}, ...
            'boxLims', config_player1{k,4}, 'relativeLimits', config_player1{k,5});
    end
    
    for k = 1:size(config_player2,1)
        if config_player2{k, 2} == 0, continue, end
        parameters_player2 = addParameter(parameters_player2, struct('model', modelCoarse_new, 'schedule', scheduleCoarse_new, 'state0', stateCoarse0_new), ...
            'name', config_player2{k,1}, 'scaling', config_player2{k,3}, ...
            'boxLims', config_player2{k,4}, 'relativeLimits', config_player2{k,5});
    end
    
else
    % Non-Nash case: Only initialize parameters for the chosen player
    if playerChoice == 1
        % Initialize Player 1 parameters
        for k = 1:size(config_player1,1)
            if config_player1{k, 2} == 0, continue, end
            parameters_player1 = addParameter(parameters_player1, struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0), ...
                'name', config_player1{k,1}, 'scaling', config_player1{k,3}, ...
                'boxLims', config_player1{k,4}, 'relativeLimits', config_player1{k,5});
        end
        % Player 2 parameters remain empty for non-Nash case
        parameters_player2 = [];
    else
        % Initialize Player 2 parameters (using config_player1 which now contains Player 2 config)
        for k = 1:size(config_player2,1)
            if config_player2{k, 2} == 0, continue, end
            parameters_player2 = addParameter(parameters_player2, struct('model', modelCoarse_new, 'schedule', scheduleCoarse_new, 'state0', stateCoarse0_new), ...
                'name', config_player2{k,1}, 'scaling', config_player2{k,3}, ...
                'boxLims', config_player2{k,4}, 'relativeLimits', config_player2{k,5});
        end
        % Player 2 parameters remain empty for non-Nash case
        parameters_player1 = [];
    end
end

%% Define mismatch functions for each twin
weighting  = {'WaterRateWeight', 1.0/(150/day), ...
              'OilRateWeight',   1.0/(80/day), ...
              'BHPWeight',      1/(20*barsa)};
% [waterWeight1, oilWeight1, bhpWeight1] = computeWeightsFromWellSol(wsRef);
% weighting = {'WaterRateWeight', [], ...
%               'OilRateWeight',   [],   ...
%               'BHPWeight',       []};
% Create objective functions that match unitBoxLM signature
player1_solver ='unitboxlm';
player2_solver = 'unitboxlm';

if NashAlgorithm ||playerChoice ==1
    if strcmp(player1_solver, 'unitboxlm')
            % Objective for Player 1 (RATE-BHP): Match RATE-BHP reference data
            mismatchFn1 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
                matchObservedOW(model, states, schedule, states_ref,...
                'computePartials', compDer, 'tstep', tstep, weighting{:},...
                'state', state, 'from_states', false, 'mismatchSum', false);
            obj_fn1 = @(p) evaluateMatchSummands(p, mismatchFn1, ...
                problem_coarse.SimulatorSetup, parameters_player1, statesRef);
    else
            mismatchFn1 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
                matchObservedOW(model, states, schedule, states_ref,...
                'computePartials', compDer, 'tstep', tstep, weighting{:},...
                'state', state, 'from_states', false);
            obj_fn1 = @(p) evaluateMatch(p, mismatchFn1, ...
                problem_coarse.SimulatorSetup, parameters_player1, statesRef);
       
    end
 end
% weighting  = {'WaterRateWeight', 0, ...
%               'OilRateWeight',   0, ...
%               'BHPWeight',      1/(20*barsa)};

if NashAlgorithm ||playerChoice ==2
    if strcmp(player2_solver, 'unitboxlm')
        % Objective for Twin 2 (BHP-RATE): Match BHP-RATE reference data
        mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
            matchObservedOW(model, states, schedule, states_ref,...
            'computePartials', compDer, 'tstep', tstep, weighting{:},...
            'state', state, 'from_states', false, 'mismatchSum', false);
        obj_fn2 = @(p) evaluateMatchSummands(p, mismatchFn2, ...
            problem_coarse_new.SimulatorSetup, parameters_player2, statesRef_new);
    else
        % Objective for Twin 2 (BHP-RATE): Match BHP-RATE reference data
        mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
            matchObservedOW(model, states, schedule, states_ref,...
            'computePartials', compDer, 'tstep', tstep, weighting{:},...
            'state', state, 'from_states', false);
        obj_fn2 = @(p) evaluateMatch(p, mismatchFn2, ...
            problem_coarse_new.SimulatorSetup, parameters_player2, statesRef_new);
    end
end

%% Create and run Nash optimization using NashOptimizationTwinProblem class

% Get initial parameter vectors
if NashAlgorithm ||playerChoice ==1
    p0_player1 = getScaledParameterVector(problem_coarse.SimulatorSetup, parameters_player1);
end

if NashAlgorithm ||playerChoice ==2
    p0_player2 = getScaledParameterVector(problem_coarse_new.SimulatorSetup, parameters_player2);
end
if NashAlgorithm
    % Create Nash optimization problem
    nashProblem = NashOptimizationTwinProblem(...
        'twinParameters', {parameters_player1, parameters_player2}, ...
        'twinProblems', {problem_coarse, problem_coarse_new}, ...
        'twinReferenceStates', {statesRef, statesRef_new}, ...
        'twinObjectives', {mismatchFn1, mismatchFn2}, ...
        'nashParameters', struct('maxIterations', 10, ...
        'convergenceTol', 1e-4, ...
        'mixingParameter', 0.9, ...
        'regularizationParam', [9.0e-6,1.5e-5], ...  % Add this
        'maxInnerIterations', 25, ...
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
    plotWellSols({wsRef, wellSols_opt1}, ...
        repmat({schedule.step.val}, 1, 2), ...
        'datasetnames', {'reference','coarse tuned RATE-BHP'});

    plotWellSols({wsRef_new, wellSols_opt2}, ...
        repmat({schedule_new.step.val}, 1, 2), ...
        'datasetnames', {'reference','coarse tuned BHP-RATE'});

    plotWellSols({wsRef_new, wellSols_opt1, wellSols_opt2}, ...
        repmat({schedule.step.val}, 1, 3), ...
        'datasetnames', {'reference','coarse tuned RATE-BHP','coarse tuned BHP-RATE'});
    %% Display final results
    fprintf('\n=== Optimization Results ===\n');
    fprintf('Final objective values:\n');
    fprintf('  RATE-BHP problem: %.4e\n', history(end).twin{1}.val(end));
    fprintf('  BHP-RATE problem: %.4e\n', history(end).twin{2}.val(end));

    fprintf('\nParameter changes from initial to final:\n');
    param_change1 = norm(u_opt{1} - p0_player1) / norm(p0_player1);
    param_change2 = norm(u_opt{2} - p0_player2) / norm(p0_player2);
    fprintf('  RATE-BHP parameters: %.2f%% change\n', param_change1 * 100);
    fprintf('  BHP-RATE parameters: %.2f%% change\n', param_change2 * 100);
else
    % Non-Nash optimization for the chosen player
    if playerChoice == 1
        fprintf('Optimizing Player 1 (RATE-BHP problem)...\n');
        if strcmp(player1_solver,'unitboxlm')
            [v1, p_opt1, history1] = unitBoxLM(p0_player1, obj_fn1, 'maxIt', 20);
        else
            [v1, p_opt1, history1] = unitBoxBFGS(p0_player1, obj_fn1, 'objChangeTol', 1e-5, ...
                                      'maxIt', 25, 'logPlot', true);
        end
        setup_opt1 = updateSetupFromScaledParameters(problem_coarse.SimulatorSetup, parameters_player1, p_opt1);
        [wellSols_opt1, states_opt1] = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, setup_opt1.schedule);
        % Plot results
        plotWellSols({wsRef, wellSols_opt1}, ...
            repmat({schedule.step.val}, 1, 2), ...
            'datasetnames', {'reference','coarse tuned RATE-BHP'});

        %% Display final results
        fprintf('\n=== Optimization Results ===\n');
        fprintf('Final objective values:\n');
        fprintf('  Player %d problem: %.4e\n', playerChoice, history1(end).val(end));
        fprintf('\nParameter changes from initial to final:\n');
        param_change1 = norm(p_opt1 - p0_player1) / norm(p0_player1);
        fprintf('  Player %d parameters: %.2f%% change\n', playerChoice, param_change1 * 100);
    else
        fprintf('Optimizing Player 2 (BHP-RATE problem)...\n');
        if strcmp(player2_solver,'unitboxlm')
            [v2, p_opt2, history2] = unitBoxLM(p0_player2, obj_fn2, 'maxIt', 20);
        else
            [v2, p_opt2, history2] = unitBoxBFGS(p0_player2, obj_fn2, 'objChangeTol', 1e-5, ...
                                      'maxIt', 25, 'logPlot', true);
        end
        setup_opt2 = updateSetupFromScaledParameters(problem_coarse_new.SimulatorSetup, parameters_player2, p_opt2);
        [wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);
        % Plot results
        plotWellSols({wsRef_new, wellSols_opt2}, ...
            repmat({schedule_new.step.val}, 1, 2), ...
            'datasetnames', {'reference','coarse tuned BHP-RATE'});

        %% Display final results
        fprintf('\n=== Optimization Results ===\n');
        fprintf('Final objective values:\n');
        fprintf('  Player %d problem: %.4e\n', playerChoice, history2(end).val(end));
        fprintf('\nParameter changes from initial to final:\n');
        param_change2 = norm(p_opt2 - p0_player2) / norm(p0_player2);
        fprintf('  Player %d parameters: %.2f%% change\n', playerChoice, param_change2 * 100);
    end

end

%% Compare reference, initial coarse and optimized coarse for a different schedule
rng(100);
W = test.schedule.control.W;

% Create validation schedule (RATE-BHP type)
s_new1 = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .2, 'perturbStep', ones(numel(dt),1));
% Simulate reference on both validation schedules
[ws_new1, states_new1] = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_new1, ...
    'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);
if NashAlgorithm || playerChoice==2
    % Create flipped validation schedule (BHP-RATE type)
    s_new2 = flipWellControlsAndAddBC(s_new1, states_new1, []);
    [ws_new2, states_new2] = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_new2, ...
        'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);
end
if NashAlgorithm
    %% Nash Algorithm Case - Validate both players
    
    % For Player 1 (RATE-BHP)
    w_opt1 = setup_opt1.schedule.control(1).W;
    w_old1 = scheduleCoarse.control(1).W;
    s_opt1 = simpleSchedule(dt, 'W', w_opt1);
    s_old1 = simpleSchedule(dt, 'W', w_old1);
    
    % For Player 2 (BHP-RATE)  
    w_opt2 = setup_opt2.schedule.control(1).W;
    w_old2 = scheduleCoarse_new.control(1).W;
    s_opt2 = simpleSchedule(dt, 'W', w_opt2);
    s_old2 = simpleSchedule(dt, 'W', w_old2);
    
    for i = 1:length(s_new2.control)
        s_opt2.control(i).W = setup_opt2.schedule.control(i).W;
        s_opt2.step.control(i) = i;
        s_old2.control(i).W = scheduleCoarse_new.control(i).W;
        s_old2.step.control(i) = i;
    end
    % Apply validation controls to both players
    for kw = 1:numel(Wref)
        % Player 1 uses RATE-BHP validation controls
        s_opt1.control.W(kw).val = s_new1.control.W(kw).val;
        s_old1.control.W(kw).val = s_new1.control.W(kw).val;
    end

    for i = 1:length(s_new2.control)
        % Player 2 uses BHP-RATE validation controls
        for kw = 1:numel(Wref)
            s_opt2.control(i).W(kw).val = s_new2.control(i).W(kw).val;
            s_old2.control(i).W(kw).val = s_new2.control(i).W(kw).val;
        end
    end
    
    % Simulate predictions for both players
    
    ws_coarse_new1 = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, s_opt1);
    ws_coarse_old1 = simulateScheduleAD(stateCoarse0, modelCoarse, s_old1);
    
    ws_coarse_new2 = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, s_opt2);
    ws_coarse_old2 = simulateScheduleAD(stateCoarse0_new, modelCoarse_new, s_old2);
    
    % Plot results for both players
    figure;
    subplot(1,2,1);
    plotWellSols({ws_new1, ws_coarse_old1, ws_coarse_new1}, ...
        repmat({schedule.step.val}, 1, 3), ...
        'datasetnames', {'reference','coarse initial','coarse tuned'});
    title('Player 1 (RATE-BHP) - Validation');
    
    subplot(1,2,2);
    plotWellSols({ws_new2, ws_coarse_old2, ws_coarse_new2}, ...
        repmat({schedule.step.val}, 1, 3), ...
        'datasetnames', {'reference','coarse initial','coarse tuned'});
    title('Player 2 (BHP-RATE) - Validation');
    
else
    %% Non-Nash Case - Validate only the chosen player
    
    if playerChoice == 1
        % Player 1 (RATE-BHP) validation
        w_opt1 = setup_opt1.schedule.control(1).W;
        w_old1 = scheduleCoarse.control(1).W;
        s_opt1 = simpleSchedule(dt, 'W', w_opt1);
        s_old1 = simpleSchedule(dt, 'W', w_old1);
        
        for kw = 1:numel(Wref)
            s_opt1.control.W(kw).val = s_new1.control.W(kw).val;
            s_old1.control.W(kw).val = s_new1.control.W(kw).val;
        end
        
        ws_coarse_new1 = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, s_opt1);
        ws_coarse_old1 = simulateScheduleAD(stateCoarse0, modelCoarse, s_old1);
        
        % Plot results for Player 1
        plotWellSols({ws_new1, ws_coarse_old1, ws_coarse_new1}, ...
            repmat({schedule.step.val}, 1, 3), ...
            'datasetnames', {'reference','coarse initial','coarse tuned'});
        title('Player 1 (RATE-BHP) - Validation');
        
    else
        % Player 2 (BHP-RATE) validation
        w_opt2 = setup_opt2.schedule.control(1).W;
        w_old2 = scheduleCoarse_new.control(1).W;
        s_opt2 = simpleSchedule(dt, 'W', w_opt2);
        s_old2 = simpleSchedule(dt, 'W', w_old2);
        

        for i = 1:length(s_new2.control)
            s_opt2.control(i).W = setup_opt2.schedule.control(i).W;
            s_opt2.step.control(i) = i;
            s_old2.control(i).W = scheduleCoarse_new.control(i).W;
            s_old2.step.control(i) = i;
        end

        for i = 1:length(s_new2.control)
            % Player 2 uses BHP-RATE validation controls
            for kw = 1:numel(Wref)
                s_opt2.control(i).W(kw).val = s_new2.control(i).W(kw).val;
                s_old2.control(i).W(kw).val = s_new2.control(i).W(kw).val;
            end
        end
        
        ws_coarse_new2 = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, s_opt2);
        ws_coarse_old2 = simulateScheduleAD(stateCoarse0_new, modelCoarse_new, s_old2);
        
        % Plot results for Player 2
        plotWellSols({ws_new2, ws_coarse_old2, ws_coarse_new2}, ...
            repmat({schedule.step.val}, 1, 3), ...
            'datasetnames', {'reference','coarse initial','coarse tuned'});
        title('Player 2 (BHP-RATE) - Validation');
    end
end
% 
% %% Store final models for later use
% if NashAlgorithm
%     final_models.nash = true;
%     final_models.player1.model = setup_opt1.model;
%     final_models.player1.state0 = setup_opt1.state0;
%     final_models.player1.schedule = setup_opt1.schedule;
%     final_models.player1.wellSols = wellSols_opt1;
%     final_models.player1.states = states_opt1;
%     
%     final_models.player2.model = setup_opt2.model;
%     final_models.player2.state0 = setup_opt2.state0;
%     final_models.player2.schedule = setup_opt2.schedule;
%     final_models.player2.wellSols = wellSols_opt2;
%     final_models.player2.states = states_opt2;
%     
%     fprintf('\n=== Final Models Stored ===\n');
%     fprintf('Both Player 1 (RATE-BHP) and Player 2 (BHP-RATE) models are available\n');
%     
% else
%     final_models.nash = false;
%     if playerChoice == 1
%         final_models.player1.model = setup_opt1.model;
%         final_models.player1.state0 = setup_opt1.state0;
%         final_models.player1.schedule = setup_opt1.schedule;
%         final_models.player1.wellSols = wellSols_opt1;
%         final_models.player1.states = states_opt1;
%         fprintf('\n=== Final Model Stored ===\n');
%         fprintf('Player 1 (RATE-BHP) model is available\n');
%     else
%         final_models.player2.model = setup_opt2.model;
%         final_models.player2.state0 = setup_opt2.state0;
%         final_models.player2.schedule = setup_opt2.schedule;
%         final_models.player2.wellSols = wellSols_opt2;
%         final_models.player2.states = states_opt2;
%         fprintf('\n=== Final Model Stored ===\n');
%         fprintf('Player 2 (BHP-RATE) model is available\n');
%     end
% end

function [waterWeight, oilWeight, bhpWeight] = computeWeightsFromWellSol(wsRef)
    % wsRef: cell array of well solutions (one entry per timestep)
    % Returns weights (1/typical value) for water rate, oil rate, BHP

    % Collect all water rates, oil rates, BHPs across all wells and all steps
    qWs_all = [];
    qOs_all = [];
    bhp_all = [];

    for step = 1:numel(wsRef)
        ws_step = wsRef{step};          % struct array for this timestep
        for w = 1:numel(ws_step)
            qWs_all = [qWs_all; ws_step(w).qWs];
            qOs_all = [qOs_all; ws_step(w).qOs];
            bhp_all = [bhp_all; ws_step(w).bhp];
        end
    end

    % Use mean absolute value as typical magnitude
    typical_qWs = mean(abs(qWs_all));
    typical_qOs = mean(abs(qOs_all));
    typical_bhp = mean(abs(bhp_all));

    % Avoid division by zero (add small epsilon if needed)
    if typical_qWs == 0, typical_qWs = 1; end
    if typical_qOs == 0, typical_qOs = 1; end
    if typical_bhp == 0, typical_bhp = 1; end

    waterWeight = 1 / typical_qWs;
    oilWeight   = 1 / typical_qOs;
    bhpWeight   = 1 / typical_bhp;
end