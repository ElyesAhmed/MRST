function [wellSols, states, scheduleReport] = ...
        simulateSequentialH2BiochemPhreeqc(state0, model, schedule, varargin)
% Run sequential-h2biochem-phreeqc with an outer same-timestep Picard iteration.
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

assert(isa(model, 'BiochemistryPhreeqcModel'), ...
    'sequential-h2biochem-phreeqc requires a BiochemistryPhreeqcModel instance.');
assert(model.phreeqcTimestepCoupling && model.isSequentialH2BiochemPhreeqcBackend(), ...
    ['simulateSequentialH2BiochemPhreeqc requires ', ...
     'phreeqcBackend=''sequential-h2biochem-phreeqc'' and phreeqcTimestepCoupling=true.']);
assert(isfield(schedule, 'step') && isfield(schedule.step, 'val') && ...
    isfield(schedule.step, 'control'), ...
    'schedule must contain step.val and step.control.');

opt = getPicardOptions(model, varargin{:});
validatePicardOptions(opt);
iterationModel = model;
iterationModel.sequentialH2BiochemPhreeqcPicardActive = true;
iterationModel = iterationModel.clearSequentialH2BiochemPhreeqcChemistryFeedback();
feedbackCleanup = onCleanup(@clearChemistryFeedback);

nSteps = numel(schedule.step.val);
assert(numel(schedule.step.control) == nSteps, ...
    'schedule.step.val and schedule.step.control must have the same length.');
[wellSols, states, reports] = deal(cell(nSteps, 1));
state = iterationModel.validateState(state0);
simulationTime = zeros(nSteps, 1);
reservoirTime = [0; cumsum(schedule.step.val(:))];

for stepNo = 1:nSteps
    printNominalStepHeader(stepNo, nSteps, reservoirTime, opt.verbose);
    timestepStart = state;
    dt = schedule.step.val(stepNo);
    validateattributes(dt, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'schedule.step.val');
    oneStepSchedule = singleStepSchedule(schedule, stepNo, nSteps);
    startReactionExtent = cumulativeReactionExtent(timestepStart, ...
        iterationModel.G.cells.num, iterationModel.biochemFluid.nbioreact);
    timer = tic();
    diagnosticState = runSequentialH2BiochemPhreeqcEquilibrium( ...
        iterationModel, timestepStart, 'inputOnly', true);
    preReactionInputCell1 = ...
        diagnosticState.sequentialH2BiochemPhreeqcPreReactionInputCell1String;

    chemistryFeedback = chemistryFeedbackFromState(iterationModel, timestepStart);
    converged = false;
    history = zeros(opt.maxIterations, 3);
    previousReactionExtent = [];
    chemistryResidual = inf;
    pHResidual = inf;
    reactionResidual = inf;
    for iteration = 1:opt.maxIterations
        if opt.showPicardProgress
            fprintf('  Picard iteration %d/%d...\n', iteration, opt.maxIterations);
        end
        % timestepStart is deliberately never overwritten in this loop.
        % The feedback is a fixed, non-AD model snapshot for this candidate
        % solve, never a nonlinear initial guess or conserved state.
        iterationModel = iterationModel.setSequentialH2BiochemPhreeqcChemistryFeedback( ...
            chemistryFeedback);
        [innerWells, innerStates, innerReport] = runInnerMrstStep( ...
            timestepStart, iterationModel, oneStepSchedule, opt);
        assert(~innerReport.Failure && numel(innerStates) == 1 && ...
            ~isempty(innerStates{1}), ...
            ['MRST failed while solving Picard iteration %d of nominal ', ...
             'timestep %d. No chemistry state was accepted.'], iteration, stepNo);

        reactedState = innerStates{1};
        reactionExtent = cumulativeReactionExtent(reactedState, ...
            iterationModel.G.cells.num, iterationModel.biochemFluid.nbioreact) - ...
            startReactionExtent;
        equilibriumState = runSequentialH2BiochemPhreeqcEquilibrium( ...
            iterationModel, reactedState);
        equilibriumState.sequentialH2BiochemPhreeqcPreReactionInputCell1String = ...
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
            if opt.showPicardProgress
                fprintf('    residuals: chemistry %.3g, pH %.3g, reaction %.3g\n', ...
                    chemistryResidual, pHResidual, reactionResidual);
            end
            if chemistryResidual <= 1 && pHResidual <= 1 && reactionResidual <= 1
                converged = true;
                break;
            end
        elseif opt.showPicardProgress
            fprintf('    initial chemistry update\n');
        end

        previousReactionExtent = reactionExtent;
        chemistryFeedback = relaxedChemistryFeedback( ...
            chemistryFeedback, equilibriumFeedback, ...
            opt.relaxation);
    end
    simulationTime(stepNo) = toc(timer);

    % if ~converged
    %     error('H2Biochem:SequentialH2BiochemPhreeqcPicardNonconvergence', ...
    %         ['sequential-h2biochem-phreeqc Picard iteration did not converge for nominal ', ...
    %          'timestep %d after %d iterations (chemistry %.3g, pH %.3g, ', ...
    %          'reaction extent %.3g; each must be <= 1).'], ...
    %         stepNo, opt.maxIterations, chemistryResidual, pHResidual, reactionResidual);
    % end

    equilibriumState.sequentialH2BiochemPhreeqcPicardIterations = iteration;
    equilibriumState.sequentialH2BiochemPhreeqcPicardChemistryResidual = chemistryResidual;
    equilibriumState.sequentialH2BiochemPhreeqcPicardPHResidual = pHResidual;
    equilibriumState.sequentialH2BiochemPhreeqcPicardReactionResidual = reactionResidual;
    equilibriumState.sequentialH2BiochemPhreeqcReactionExtentMoles = reactionExtent;
    equilibriumState.sequentialH2BiochemPhreeqcPicardHistory = history(2:iteration, :);

    state = equilibriumState;
    states{stepNo} = state;
    wellSols{stepNo} = innerWells{1};
    controlReport = innerReport.ControlstepReports{1};
    controlReport.SequentialH2BiochemPhreeqcPicard = struct( ...
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

fprintf('*** Simulation complete. Solved %d control steps in %s ***\n', ...
    nSteps, formatTimeRange(sum(simulationTime)));

clearChemistryFeedback();
clear feedbackCleanup

    function clearChemistryFeedback()
        iterationModel = iterationModel.clearSequentialH2BiochemPhreeqcChemistryFeedback();
        iterationModel.sequentialH2BiochemPhreeqcPicardActive = false;
    end
end

function opt = getPicardOptions(model, varargin)
configured = model.phreeqcCouplingOptions;
opt = struct( ...
    'maxIterations', configuredValue(configured, 'sequentialH2BiochemPhreeqcMaxIterations', 30), ...
    'relaxation', configuredValue(configured, 'sequentialH2BiochemPhreeqcRelaxation', 0.25), ...
    'absoluteTolerance', configuredValue(configured, ...
        'sequentialH2BiochemPhreeqcAbsoluteTolerance', 1e-8), ...
    'relativeTolerance', configuredValue(configured, ...
        'sequentialH2BiochemPhreeqcRelativeTolerance', 1e-2), ...
    'pHTolerance', configuredValue(configured, 'sequentialH2BiochemPhreeqcPHTolerance', 2e-2), ...
    'pKaTolerance', configuredValue(configured, 'sequentialH2BiochemPhreeqcPkaTolerance', 2e-2), ...
    'reactionAbsoluteTolerance', configuredValue(configured, ...
        'sequentialH2BiochemPhreeqcReactionAbsoluteTolerance', 1e-12), ...
    'reactionRelativeTolerance', configuredValue(configured, ...
        'sequentialH2BiochemPhreeqcReactionRelativeTolerance', 1e-3), ...
    'nonlinearSolver', [], ...
    'linearSolver', [], ...
    'verbose', mrstVerbose(), ...
    'showPicardProgress', true, ...
    'suppressInnerOutput', true, ...
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
validateattributes(opt.showPicardProgress, {'logical', 'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, 'showPicardProgress');
validateattributes(opt.suppressInnerOutput, {'logical', 'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, 'suppressInnerOutput');
end

function [wellSols, states, report] = runInnerMrstStep( ...
        state0, model, schedule, opt)
    args = { ...
        'Verbose', opt.verbose, ...
        'OutputMinisteps', false, ...
        'NonLinearSolver', opt.nonlinearSolver, ...
        'LinearSolver', opt.linearSolver, ...
        'checkOperators', opt.checkOperators};
    if opt.suppressInnerOutput
        [capturedOutput, wellSols, states, report] = evalc( ... %#ok<ASGLU>
            'simulateScheduleAD(state0, model, schedule, args{:})');
        capturedOutput = regexprep(capturedOutput, ...
            'Solving timestep[^\r\n]*\r?\n?', '');
        capturedOutput = regexprep(capturedOutput, ...
            '\*\*\* Simulation complete\.[^\r\n]*\r?\n?', '');
        if ~isempty(strtrim(capturedOutput))
            fprintf('%s', capturedOutput);
        end
    else
        [wellSols, states, report] = simulateScheduleAD( ...
            state0, model, schedule, args{:});
    end
end

function printNominalStepHeader(stepNo, nSteps, reservoirTime, verbose)
    nDigits = floor(log10(nSteps)) + 1;
    if verbose
        nChar = 0;
    else
        nChar = numel(formatTimeRange(reservoirTime(end), 2));
    end
    fprintf('Solving timestep %0*d/%0*d: %-*s -> %s\n', ...
        nDigits, stepNo, nDigits, nSteps, nChar, ...
        formatTimeRange(reservoirTime(stepNo), 2), ...
        formatTimeRange(reservoirTime(stepNo + 1), 2));
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
assert(~isempty(index), 'sequential-h2biochem-phreeqc requires an EOS CO2 component.');
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
    'sequential-h2biochem-phreeqc feedback field %s must be a finite cell vector.', field);
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
if isfield(state, 'sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles')
    extent = value(state.sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles);
else
    extent = zeros(nc, nreact);
end
assert(isequal(size(extent), [nc, nreact]), ...
    'sequential-h2biochem-phreeqc cumulative H2 consumption has invalid dimensions.');
assert(all(isfinite(extent(:))), ...
    'h2-biochem cumulative reaction extent contains non-finite values.');
end
