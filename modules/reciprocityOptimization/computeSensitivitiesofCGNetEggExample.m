%% Automated Parameter Partitioning using Sensitivity Analysis
mrstModule add ad-core ad-blackoil deckformat agglom upscaling coarsegrid...
        mrst-gui ad-props incomp optimization test-suite linearsolvers

%% Setup reference model (same as before)
test    = TestCase('egg_wo', 'realization', 1);
problem = test.getPackedSimulationProblem();
problem.SimulatorSetup.model.useCNVConvergence = false;
problem.SimulatorSetup.model.nonlinearTolerance = 1.0e-9;
problem.SimulatorSetup.model.toleranceMB = 1.0e-9;
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

% Create schedule with flipped well controls
schedule_new = flipWellControlsAndAddBC(schedule, statesRef, []);
problem_new  = packSimulationProblem(problem.SimulatorSetup.state0, modelRef, schedule_new, 'egg_wo_new');
problem_new.SimulatorSetup.model = modelRef;
simulatePackedProblem(problem_new);
[wsRef_new, statesRef_new] = getPackedSimulatorOutput(problem_new);

%% Create coarse models for both problems
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

% Coarse model for BHP-Inj-RATE-Prod problem  
modelCoarse_new = upscaleModelTPFA(modelRef, blockIx);
modelCoarse_new.AutoDiffBackend = AutoDiffBackend();
pts_new = modelCoarse_new.fluid.krPts;
scaling_new = {'SWL', pts_new.w(1), 'SWCR', pts_new.w(2), 'SWU', pts_new.w(3), ...
           'SOWCR', pts_new.ow(2), 'KRW', pts_new.w(4), 'KRO', pts_new.ow(4)};
modelCoarse_new = imposeRelpermScaling(modelCoarse_new, scaling_new{:});
modelCoarse_new.toleranceCNV = 1e-6;

%% Set up initial states and schedules
stateCoarse0 = upscaleState(modelCoarse, modelRef, test.state0);
scheduleCoarse = upscaleSchedule(modelCoarse, schedule, 'wellUpscaleMethod', 'sum');

stateCoarse0_new = upscaleState(modelCoarse_new, modelRef, test.state0);
scheduleCoarse_new = upscaleSchedule(modelCoarse_new, schedule_new, 'wellUpscaleMethod', 'sum');

%% Simulate initial coarse models
[wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);
[wsCoarse_new, statesCoarse_new] = simulateScheduleAD(stateCoarse0_new, modelCoarse_new, scheduleCoarse_new);

%% Define ALL parameters for sensitivity analysis (include all for both problems)
pv = modelCoarse.operators.pv;
pv_new = modelCoarse_new.operators.pv;

% Common parameter configuration (ALL parameters active)
config_all = {...
    %name           include    scaling              boxlims   relativeLimits
    'porevolume',       1,     'linear',    [.005*pv, 2.5*pv],              []
    'conntrans',        1,        'log',                  [],     [1e-3, 1e3]
    'transmissibility', 1,        'log',                  [],     [1e-3, 1e3]
    'swl',              1,     'linear',             [0, .5],              []
    'swcr',             1,     'linear',             [0, .5],              []
    'swu',              1,     'linear',             [.5, 1],              []
    'sowcr',            1,     'linear',             [0, .5],              []
    'krw',              1,     'linear',           [.25, 2.5],              []
    'kro',              1,     'linear',           [.25, 2.5],              []
    'sw',               1,     'linear',            [.005, 0.5],            []
    'pressure',         1,     'log',           [.895.*stateCoarse0.pressure, 1.5.*stateCoarse0.pressure], []
    };

%% Initialize ALL parameters for both problems
parameters_all1 = [];
for k = 1:size(config_all,1)
    parameters_all1 = addParameter(parameters_all1, struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0), ...
        'name', config_all{k,1}, 'scaling', config_all{k,3}, ...
        'boxLims', config_all{k,4}, 'relativeLimits', config_all{k,5});
end

parameters_all2 = [];
for k = 1:size(config_all,1)
    parameters_all2 = addParameter(parameters_all2, struct('model', modelCoarse_new, 'schedule', scheduleCoarse_new, 'state0', stateCoarse0_new), ...
        'name', config_all{k,1}, 'scaling', config_all{k,3}, ...
        'boxLims', config_all{k,4}, 'relativeLimits', config_all{k,5});
end

%% Define mismatch functions for sensitivity analysis
weighting  = {'WaterRateWeight', 1.0/(150/day), ...
              'OilRateWeight',   1.0/(80/day), ...
              'BHPWeight',      1/(20*barsa)};

% Objective function handles for sensitivity analysis
mismatchFn1_sens = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
    'computePartials', compDer, 'tstep', tstep, weighting{:},...
    'state', state, 'from_states', false);

mismatchFn2_sens = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
    'computePartials', compDer, 'tstep', tstep, weighting{:},...
    'state', state, 'from_states', false);

%% Compute sensitivities for Problem 1 (RATE-BHP)
fprintf('Computing sensitivities for RATE-BHP problem...\n');
setup1 = struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0);
objh1 = @(tstep, model, state) mismatchFn1_sens(model, statesCoarse, scheduleCoarse, statesRef,...
    true, tstep, state);

raw_sens1 = computeSensitivitiesAdjointAD(setup1, statesCoarse, parameters_all1, objh1);

% Process and scale sensitivities for Problem 1
sens1 = struct();
for kp = 1:numel(parameters_all1)
    nm = parameters_all1{kp}.name;
    if strcmp(parameters_all1{kp}.type, 'multiplier')
        sens1.(nm) = (raw_sens1.(nm)) .* parameters_all1{kp}.getParameterValue(setup1);
    else
        sens1.(nm) = raw_sens1.(nm);
    end
end

%% Compute sensitivities for Problem 2 (BHP-RATE)
fprintf('Computing sensitivities for BHP-RATE problem...\n');
setup2 = struct('model', modelCoarse_new, 'schedule', scheduleCoarse_new, 'state0', stateCoarse0_new);
objh2 = @(tstep, model, state) mismatchFn2_sens(model, statesCoarse_new, scheduleCoarse_new, statesRef_new,...
    true, tstep, state);

raw_sens2 = computeSensitivitiesAdjointAD(setup2, statesCoarse_new, parameters_all2, objh2);

% Process and scale sensitivities for Problem 2
sens2 = struct();
for kp = 1:numel(parameters_all2)
    nm = parameters_all2{kp}.name;
    if strcmp(parameters_all2{kp}.type, 'multiplier')
        sens2.(nm) = (raw_sens2.(nm)) .* parameters_all2{kp}.getParameterValue(setup2);
    else
        sens2.(nm) = raw_sens2.(nm);
    end
end

%% Analyze and compare sensitivities to determine optimal partitioning
fprintf('\n=== SENSITIVITY ANALYSIS FOR PARAMETER PARTITIONING ===\n');

% Extract parameter names from ModelParameter objects
param_names = cell(1, numel(parameters_all1));
for k = 1:numel(parameters_all1)
    param_names{k} = parameters_all1{k}.name;
end

% Extract sensitivity magnitudes (L2 norm)
sensitivity_magnitudes = zeros(numel(param_names), 2);

for k = 1:numel(param_names)
    nm = param_names{k};
    % Check if sensitivity field exists and extract magnitude
    if isfield(sens1, nm) && ~isempty(sens1.(nm))
        sensitivity_magnitudes(k, 1) = norm(sens1.(nm));  % Problem 1 sensitivity
    else
        sensitivity_magnitudes(k, 1) = 0;
    end
    
    if isfield(sens2, nm) && ~isempty(sens2.(nm))
        sensitivity_magnitudes(k, 2) = norm(sens2.(nm));  % Problem 2 sensitivity
    else
        sensitivity_magnitudes(k, 2) = 0;
    end
end

% Calculate sensitivity ratios (avoid division by zero)
sensitivity_ratios = zeros(size(sensitivity_magnitudes, 1), 1);
for k = 1:size(sensitivity_magnitudes, 1)
    if sensitivity_magnitudes(k, 2) > 0
        sensitivity_ratios(k) = sensitivity_magnitudes(k, 1) / sensitivity_magnitudes(k, 2);
    else
        sensitivity_ratios(k) = inf; % Infinite ratio if P2 sensitivity is zero
    end
end

% Determine optimal partitioning based on sensitivity dominance
player1_params = {};
player2_params = {};
shared_params = {};

sensitivity_threshold = 2.0; % Parameters with ratio > 2 go to Player 1, < 0.5 go to Player 2

fprintf('\nParameter Sensitivity Analysis:\n');
fprintf('%-20s %-15s %-15s %-10s %-10s\n', 'Parameter', 'Sens_P1', 'Sens_P2', 'Ratio', 'Assignment');
fprintf([repmat('-', 1, 70) '\n']);

for k = 1:numel(param_names)
    nm = param_names{k};
    sens_p1 = sensitivity_magnitudes(k, 1);
    sens_p2 = sensitivity_magnitudes(k, 2);
    ratio = sensitivity_ratios(k);
    
    if ratio > sensitivity_threshold
        assignment = 'Player 1';
        player1_params{end+1} = nm;
    elseif ratio < 1/sensitivity_threshold && ~isinf(ratio)
        assignment = 'Player 2';
        player2_params{end+1} = nm;
    else
        assignment = 'Shared';
        shared_params{end+1} = nm;
    end
    
    if isinf(ratio)
        ratio_str = 'Inf';
    else
        ratio_str = sprintf('%.2f', ratio);
    end
    
    fprintf('%-20s %-15.2e %-15.2e %-10s %-10s\n', nm, sens_p1, sens_p2, ratio_str, assignment);
end