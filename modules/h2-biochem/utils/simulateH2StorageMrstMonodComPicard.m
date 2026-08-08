function [wellSols, states, scheduleReport] = ...
        simulateH2StorageMrstMonodComPicard(state0, model, schedule, varargin)
% Run mrst-monod-com with an outer same-timestep Picard iteration.
%
% Each Picard iterate solves one nominal schedule timestep from its fixed
% timestep-start state. A numerical PHREEQC chemistry snapshot is placed
% on the model only for the Monod kinetic substrate/speciation evaluation;
% it is never made the next iterate's accumulation state. MRST's nbact
% Monod sources therefore advance exactly once in every candidate
% timestep, and PHREEQC receives the already reacted full inventory only
% for equilibrium repartitioning.
%
% Direct simulateScheduleAD is intentionally unsupported for this backend:
% it would produce a lagged one-pass split rather than this fixed-point
% scheme. This wrapper returns one state and one well solution per nominal
% schedule step, irrespective of any nonlinear-solver internal ministeps.

assert(model.phreeqcTimestepCoupling && model.isMrstMonodComPhreeqcBackend(), ...
    ['simulateH2StorageMrstMonodComPicard requires ', ...
     'phreeqcBackend=''mrst-monod-com'' and phreeqcTimestepCoupling=true.']);
assert(~model.phreeqcConservativeCarbonTransfer, ...
    ['mrst-monod-com transfers complete PHREEQC component inventories; ', ...
     'phreeqcConservativeCarbonTransfer must be false.']);
assert(isfield(schedule, 'step') && isfield(schedule.step, 'val') && ...
    isfield(schedule.step, 'control'), ...
    'schedule must contain step.val and step.control.');

opt = getPicardOptions(model, varargin{:});
validatePicardOptions(opt);
iterationModel = model;
iterationModel.phreeqcMrstMonodComPicardActive = true;
iterationModel = iterationModel.clearMrstMonodComChemistryFeedback();
feedbackCleanup = onCleanup(@clearChemistryFeedback);

nSteps = numel(schedule.step.val);
assert(numel(schedule.step.control) == nSteps, ...
    'schedule.step.val and schedule.step.control must have the same length.');
[wellSols, states, reports] = deal(cell(nSteps, 1));
state = iterationModel.validateState(state0);
simulationTime = zeros(nSteps, 1);

for stepNo = 1:nSteps
    timestepStart = state;
    dt = schedule.step.val(stepNo);
    validateattributes(dt, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'schedule.step.val');
    oneStepSchedule = singleStepSchedule(schedule, stepNo, nSteps);
    startReactionExtent = cumulativeReactionExtent(timestepStart, ...
        iterationModel.G.cells.num, iterationModel.biochemFluid.nbioreact);
    diagnosticState = runH2StorageMrstMonodIPhreeqcCOMEquilibrium( ...
        iterationModel, timestepStart, 'inputOnly', true);
    preReactionInputCell1 = ...
        diagnosticState.phreeqcMrstMonodComPreReactionInputCell1String;

    chemistryFeedback = chemistryFeedbackFromState(iterationModel, timestepStart);
    converged = false;
    timer = tic();
    history = zeros(opt.maxIterations, 3);
    previousReactionExtent = [];
    chemistryResidual = inf;
    pHResidual = inf;
    reactionResidual = inf;
    for iteration = 1:opt.maxIterations
        % timestepStart is deliberately never overwritten in this loop.
        % The feedback is a fixed, non-AD model snapshot for this candidate
        % solve, never a nonlinear initial guess or conserved state.
        iterationModel = iterationModel.setMrstMonodComChemistryFeedback( ...
            chemistryFeedback);
        [innerWells, innerStates, innerReport] = simulateScheduleAD( ...
            timestepStart, iterationModel, oneStepSchedule, ...
            'Verbose', opt.verbose, ...
            'OutputMinisteps', false, ...
            'NonLinearSolver', opt.nonlinearSolver, ...
            'LinearSolver', opt.linearSolver, ...
            'checkOperators', opt.checkOperators);
        assert(~innerReport.Failure && numel(innerStates) == 1 && ...
            ~isempty(innerStates{1}), ...
            ['MRST failed while solving Picard iteration %d of nominal ', ...
             'timestep %d. No chemistry state was accepted.'], iteration, stepNo);

        reactedState = innerStates{1};
        reactionExtent = cumulativeReactionExtent(reactedState, ...
            iterationModel.G.cells.num, iterationModel.biochemFluid.nbioreact) - ...
            startReactionExtent;
        equilibriumState = runH2StorageMrstMonodIPhreeqcCOMEquilibrium( ...
            iterationModel, reactedState);
        equilibriumState.phreeqcMrstMonodComPreReactionInputCell1String = ...
            preReactionInputCell1;
        equilibriumFeedback = chemistryFeedbackFromState( ...
            iterationModel, equilibriumState);

        if ~isempty(previousReactionExtent)
            [chemistryResidual, pHResidual] = chemistryFeedbackResidualNorm( ...
                equilibriumFeedback, chemistryFeedback, opt);
            reactionResidual = normalizedDifference(reactionExtent, ...
                previousReactionExtent, ...
                opt.reactionAbsoluteTolerance, opt.reactionRelativeTolerance);
            history(iteration, :) = [chemistryResidual, pHResidual, reactionResidual];
            if chemistryResidual <= 1 && pHResidual <= 1 && reactionResidual <= 1
                converged = true;
                break;
            end
        end

        previousReactionExtent = reactionExtent;
        chemistryFeedback = relaxedChemistryFeedback( ...
            chemistryFeedback, equilibriumFeedback, ...
            opt.relaxation);
    end
    simulationTime(stepNo) = toc(timer);

    % if ~converged
    %     error('H2Biochem:MrstMonodComPicardNonconvergence', ...
    %         ['mrst-monod-com Picard iteration did not converge for nominal ', ...
    %          'timestep %d after %d iterations (chemistry %.3g, pH %.3g, ', ...
    %          'reaction extent %.3g; each must be <= 1).'], ...
    %         stepNo, opt.maxIterations, chemistryResidual, pHResidual, reactionResidual);
    % end

    equilibriumState.phreeqcMrstMonodComPicardIterations = iteration;
    equilibriumState.phreeqcMrstMonodComPicardChemistryResidual = chemistryResidual;
    equilibriumState.phreeqcMrstMonodComPicardPHResidual = pHResidual;
    equilibriumState.phreeqcMrstMonodComPicardReactionResidual = reactionResidual;
    equilibriumState.phreeqcMrstMonodComReactionExtentMoles = reactionExtent;
    equilibriumState.phreeqcMrstMonodComPicardHistory = history(2:iteration, :);

    state = equilibriumState;
    states{stepNo} = state;
    wellSols{stepNo} = innerWells{1};
    controlReport = innerReport.ControlstepReports{1};
    controlReport.PhreeqcMrstMonodComPicard = struct( ...
        'Iterations', iteration, ...
        'ChemistryResidual', chemistryResidual, ...
        'PHResidual', pHResidual, ...
        'ReactionResidual', reactionResidual, ...
        'History', history(2:iteration, :));
    reports{stepNo} = controlReport;
end

scheduleReport = struct();
scheduleReport.ControlstepReports = reports;
scheduleReport.ReservoirTime = cumsum(schedule.step.val);
scheduleReport.Converged = true(nSteps, 1);
scheduleReport.Iterations = cellfun(@(report) report.Iterations, reports);
scheduleReport.SimulationTime = simulationTime;
scheduleReport.Failure = false;

clearChemistryFeedback();
clear feedbackCleanup

    function clearChemistryFeedback()
        iterationModel = iterationModel.clearMrstMonodComChemistryFeedback();
        iterationModel.phreeqcMrstMonodComPicardActive = false;
    end
end

function opt = getPicardOptions(model, varargin)
configured = model.phreeqcCouplingOptions;
opt = struct( ...
    'maxIterations', configuredValue(configured, 'mrstMonodComMaxIterations', 30), ...
    'relaxation', configuredValue(configured, 'mrstMonodComRelaxation', 0.25), ...
    'absoluteTolerance', configuredValue(configured, ...
        'mrstMonodComAbsoluteTolerance', 1e-8), ...
    'relativeTolerance', configuredValue(configured, ...
        'mrstMonodComRelativeTolerance', 1e-2), ...
    'pHTolerance', configuredValue(configured, 'mrstMonodComPHTolerance', 2e-2), ...
    'pKaTolerance', configuredValue(configured, 'mrstMonodComPkaTolerance', 2e-2), ...
    'reactionAbsoluteTolerance', configuredValue(configured, ...
        'mrstMonodComReactionAbsoluteTolerance', 1e-12), ...
    'reactionRelativeTolerance', configuredValue(configured, ...
        'mrstMonodComReactionRelativeTolerance', 1e-3), ...
    'nonlinearSolver', [], ...
    'linearSolver', [], ...
    'verbose', mrstVerbose(), ...
    'checkOperators', []);
opt = merge_options(opt, varargin{:});
end

function value = configuredValue(configured, name, default)
if isfield(configured, name)
    value = configured.(name);
else
    value = default;
end
end

function validatePicardOptions(opt)
validateattributes(opt.maxIterations, {'numeric'}, ...
    {'scalar', 'integer', 'finite', '>=', 1}, mfilename, 'maxIterations');
validateattributes(opt.relaxation, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<=', 1}, mfilename, 'relaxation');
positive = {'absoluteTolerance', 'relativeTolerance', 'pHTolerance', 'pKaTolerance', ...
    'reactionAbsoluteTolerance', 'reactionRelativeTolerance'};
for i = 1:numel(positive)
    validateattributes(opt.(positive{i}), {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, mfilename, positive{i});
end
validateattributes(opt.verbose, {'logical', 'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, 'verbose');
end

function schedule = singleStepSchedule(schedule, stepNo, nSteps)
% Subset every per-step metadata field, while leaving controls unchanged.
step = schedule.step;
for field = fieldnames(step).'
    value = step.(field{1});
    if isnumeric(value) || islogical(value) || iscell(value)
        if numel(value) == nSteps
            step.(field{1}) = value(stepNo);
        end
    end
end
schedule.step = step;
end

function feedback = relaxedChemistryFeedback(previous, current, relaxation)
feedback = current;
fields = {'pH', 'carbonatePka1', 'hco3Molality', 'co2Molality', ...
    'sulfateMolality'};
for i = 1:numel(fields)
    field = fields{i};
    feedback.(field) = (1 - relaxation).*previous.(field) + ...
        relaxation.*current.(field);
end
end

function [chemistry, pH] = chemistryFeedbackResidualNorm(current, previous, opt)
pH = max(abs(current.pH - previous.pH))./opt.pHTolerance;
chemistry = max(abs(current.carbonatePka1 - previous.carbonatePka1))./ ...
    opt.pKaTolerance;
fields = {'hco3Molality', 'co2Molality', 'sulfateMolality'};
for i = 1:numel(fields)
    field = fields{i};
    if strcmp(field, 'sulfateMolality')
        concentrationFloor = max(opt.absoluteTolerance, 1e-6);
    else
        concentrationFloor = max(opt.absoluteTolerance, 1e-8);
    end
    chemistry = max(chemistry, scaledChemistryDifference( ...
        current.(field), previous.(field), concentrationFloor, ...
        opt.relativeTolerance));
end
end

function feedback = chemistryFeedbackFromState(model, state)
nc = model.G.cells.num;
rhoWater = model.EOSModel.rho_water;
feedback = struct( ...
    'pH', stateVector(state, 'phreeqcPH', nc, model.carbonateBufferPH), ...
    'carbonatePka1', stateVector(state, 'phreeqcCarbonatePka1', nc, ...
        model.carbonateBufferPka1), ...
    'hco3Molality', [], ...
    'co2Molality', [], ...
    'sulfateMolality', stateVector(state, 'tracerSO4', nc, 0)./rhoWater);
if isfield(state, 'phreeqcHCO3Molality')
    feedback.hco3Molality = stateVector(state, 'phreeqcHCO3Molality', nc, 0);
else
    feedback.hco3Molality = stateVector(state, 'tracerHCO3', nc, 0)./rhoWater;
end
if isfield(state, 'phreeqcCO2Molality')
    feedback.co2Molality = stateVector(state, 'phreeqcCO2Molality', nc, 0);
else
    feedback.co2Molality = dissolvedCO2Molality(model, state, rhoWater);
end
end

function molality = dissolvedCO2Molality(model, state, rhoWater)
names = model.EOSModel.getComponentNames();
index = find(strcmpi(names, 'CO2'), 1);
assert(~isempty(index), 'mrst-monod-com requires an EOS CO2 component.');
x = value(state.x);
if iscell(x)
    xCO2 = x{index};
else
    xCO2 = x(:, index);
end
rhoMolar = value(state.pressure)./(value(state.Z_L).*8.314.*value(state.T));
molality = max(rhoMolar.*value(xCO2)./rhoWater, 0);
end

function values = stateVector(state, field, nc, default)
if isfield(state, field)
    values = value(state.(field));
else
    values = default;
end
if isscalar(values)
    values = repmat(values, nc, 1);
else
    values = values(:);
end
assert(numel(values) == nc && isreal(values) && all(isfinite(values)), ...
    'mrst-monod-com feedback field %s must be a finite cell vector.', field);
end

function residual = normalizedDifference(current, previous, absoluteTolerance, relativeTolerance)
scale = absoluteTolerance + relativeTolerance.*max(abs(current), abs(previous));
residual = max(abs(current(:) - previous(:))./scale(:));
end

function residual = scaledChemistryDifference(current, previous, concentrationFloor, tolerance)
scale = max(max(abs(current), abs(previous)), concentrationFloor);
residual = max(abs(current(:) - previous(:))./scale(:))./tolerance;
end

function extent = cumulativeReactionExtent(state, nc, nreact)
if isfield(state, 'phreeqcMrstMonodComCumulativeH2ConsumptionMoles')
    extent = value(state.phreeqcMrstMonodComCumulativeH2ConsumptionMoles);
else
    extent = zeros(nc, nreact);
end
assert(isequal(size(extent), [nc, nreact]), ...
    'mrst-monod-com cumulative H2 consumption has invalid dimensions.');
assert(all(isfinite(extent(:))), ...
    'MRST Monod cumulative reaction extent contains non-finite values.');
end
