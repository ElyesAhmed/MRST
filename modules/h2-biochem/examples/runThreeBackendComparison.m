function scenarios = runThreeBackendComparison(varargin)
% Compare original, compositional-PHREEQC, and h2-biochem-PHREEQC backends.
%
% Aqueous H2 and CO2 profiles are shown at the end of injection, storage,
% and production. Additional figures compare H2 loss over time and its final
% spatial distribution. Results remain in memory and no files are written.
%
% EXAMPLE:
%   db = '\\wsl.localhost\Ubuntu\path\to\PHREEQC_Modified.DAT';
%   scenarios = runThreeBackendComparison('phreeqcDatabaseFile', db);
%   scenarios = runThreeBackendComparison('phreeqcDatabaseFile', db, ...
%       'referenceUseSoreideWhitsonEOS', true);
% The hybrid case enables adaptive Picard relaxation and cuts only the
% schedule step that fails the outer coupling convergence test.

    mrstModule add ad-props compositional deckformat h2-biochem

    opt = struct( ...
        'phreeqcDatabaseFile', '', ...
        'phreeqcPicardMaxIterations', 15, ...
        'phreeqcPicardRelaxation', 0.75, ...
        'hybridAdaptiveRelaxation', true, ...
        'hybridFailOnNonconvergence', true, ...
        'hybridCutTimestepOnNonconvergence', true, ...
        'hybridMinimumRelaxation', 0.4, ...
        'hybridMaxTimestepCuts', 2, ...
        'referenceUseSoreideWhitsonEOS', false);
    opt = merge_options(opt, varargin{:});
    validateattributes(opt.referenceUseSoreideWhitsonEOS, ...
        {'logical', 'numeric'}, {'scalar', 'real', 'finite'}, ...
        mfilename, 'referenceUseSoreideWhitsonEOS');
    opt.referenceUseSoreideWhitsonEOS = ...
        logical(opt.referenceUseSoreideWhitsonEOS);
    validateattributes(opt.hybridMinimumRelaxation, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive', ...
        '<=', opt.phreeqcPicardRelaxation}, ...
        mfilename, 'hybridMinimumRelaxation');
    validateattributes(opt.hybridMaxTimestepCuts, {'numeric'}, ...
        {'scalar', 'integer', 'finite', 'nonnegative'}, ...
        mfilename, 'hybridMaxTimestepCuts');
    databaseFile = resolvePhreeqcDatabaseFile(opt.phreeqcDatabaseFile);

    commonOptions = { ...
        'rate', 'highrate', ...
        'scheduleMode', 'complete', ...
        'injectionCO2', 0, ...
        'bacteriamodel', true, ...
        'bactDiffusion', false, ...
        'chemotaxisEffect', false, ...
        'molecularDiffusion', false, ...
        'molecularDispersion', false, ...
        'bioClogging', false, ...
        'carbonateBuffer', true, ...
        'carbonateBufferPH', 6.24, ...
        'initialHCO3', 1.370e-3, ...
        'initialOverallCO2', 0};

    cases = [ ...
        makeCase("Original h2-biochem", { ...
            'equilibrateInitialCO2', false, ...
            'phreeqcTimestepCoupling', false}, false), ...
        makeCase("Sequential h2-biochem PHREEQC", { ...
            'equilibrateInitialCO2', false, ...
            'paperBiomassKinetics', true, ...
            'nbact0', 1, ...
            'phreeqcTimestepCoupling', true, ...
            'phreeqcBackend', 'sequential-h2biochem-phreeqc', ...
            'phreeqcPicardRelativeTolerance', 5e-2, ...
            'phreeqcPicardReactionRelativeTolerance', 5e-2, ...
            'phreeqcPicardPHTolerance', 2e-2, ...
            'phreeqcPicardPkaTolerance', 2e-2, ...
            'phreeqcDatabaseFile', databaseFile, ...
            'phreeqcPicardMaxIterations', opt.phreeqcPicardMaxIterations, ...
            'phreeqcPicardRelaxation', opt.phreeqcPicardRelaxation}, true), ...
        makeCase("Sequential compositional PHREEQC", { ...
            'equilibrateInitialCO2', false, ...
            'phreeqcTimestepCoupling', true, ...
            'phreeqcBackend', 'sequential-compositional-phreeqc', ...
            'phreeqcDatabaseFile', databaseFile}, false)];

    scenarios = initializeScenarios(numel(cases));
    for caseNo = 1:numel(cases)
        current = cases(caseNo);
        fprintf('\nBackend %d/%d: %s\n', ...
            caseNo, numel(cases), current.name);
        [~, model, schedule, state0] = ...
            setupH2StorageExampleWithSRB_benchmark( ...
            commonOptions{:}, current.options{:});
        if current.sequentialDriver
            assert(model.bacterialDecayOrder == 1, ...
               'The UGFACT-matched hybrid case requires first-order biomass decay.');
            assert(max(abs(value(state0.nbact) - 1), [], 'all') < 1e-12, ...
               'The UGFACT-matched hybrid case requires initial N/N0 = 1.');
        end

        solver = NonLinearSolver();
        solver.maxTimestepCuts = 12;
        timer = tic();
        if current.sequentialDriver
            [ws, states, report] = simulateSequentialH2BiochemPhreeqc( ...
                state0, model, schedule, ...
                'nonlinearSolver', solver, ...
                'adaptiveRelaxation', opt.hybridAdaptiveRelaxation, ...
                'minimumRelaxation', opt.hybridMinimumRelaxation, ...
                'failOnNonconvergence', opt.hybridFailOnNonconvergence, ...
                'cutTimestepOnNonconvergence', ...
                    opt.hybridCutTimestepOnNonconvergence, ...
                'maxPicardTimestepCuts', opt.hybridMaxTimestepCuts);
        else
            [ws, states, report] = simulateScheduleAD( ...
                state0, model, schedule, 'nonlinearSolver', solver);
        end
        if isstruct(report) && isfield(report, 'Failure')
            assert(~report.Failure, ...
                'Simulation report indicates failure for "%s".', current.name);
        end

        metrics = collectMetrics(states, schedule, model);
        scenarios(caseNo).name = current.name;
        scenarios(caseNo).model = model;
        scenarios(caseNo).schedule = schedule;
        scenarios(caseNo).states = states;
        scenarios(caseNo).ws = ws;
        scenarios(caseNo).state0 = state0;
        scenarios(caseNo).timeDays = metrics.timeDays;
        scenarios(caseNo).snapshotIndices = metrics.snapshotIndices;
        scenarios(caseNo).snapshotLabels = metrics.snapshotLabels;
        scenarios(caseNo).aqueousH2 = metrics.aqueousH2;
        scenarios(caseNo).aqueousCO2 = metrics.aqueousCO2;
        scenarios(caseNo).lossPercent = metrics.lossPercent;
        scenarios(caseNo).finalSpatialConsumptionMoles = ...
            metrics.finalSpatialConsumptionMoles;
        scenarios(caseNo).xDimensionless = metrics.xDimensionless;
        scenarios(caseNo).injectedH2Moles = metrics.injectedH2Moles;
        scenarios(caseNo).consumedH2Moles = metrics.consumedH2Moles;
        scenarios(caseNo).runtimeSeconds = toc(timer);
    end

    ugfactRoot = fileparts(fileparts(fileparts(databaseFile)));
    scenarios(end + 1) = runUGFACTReference( ...
        ugfactRoot, opt.referenceUseSoreideWhitsonEOS);
    printSummary(scenarios);
    createComparisonFigures(scenarios);
end

function item = makeCase(name, options, sequentialDriver)
    item = struct('name', name, 'options', {options}, ...
        'sequentialDriver', sequentialDriver);
end

function scenarios = initializeScenarios(nCases)
    template = struct( ...
        'name', "", 'model', [], 'schedule', [], 'states', [], ...
        'ws', [], 'state0', [], ...
        'timeDays', [], 'snapshotIndices', [], 'snapshotLabels', [], ...
        'aqueousH2', [], 'aqueousCO2', [], 'lossPercent', [], ...
        'finalSpatialConsumptionMoles', [], 'xDimensionless', [], ...
        'injectedH2Moles', nan, 'consumedH2Moles', nan, ...
        'runtimeSeconds', nan);
    scenarios = repmat(template, nCases, 1);
end

function metrics = collectMetrics(states, schedule, model)
    assert(numel(states) == numel(schedule.step.val), ...
        'Saved states and schedule must contain the same number of steps.');
    componentNames = model.EOSModel.getComponentNames();
    idxH2 = findComponentIndex(componentNames, {'H2', 'Hydrogen'});
    idxCO2 = findComponentIndex(componentNames, {'CO2', 'CarbonDioxide'});

    controls = schedule.step.control(:);
    snapshotIndices = [find(controls == 1, 1, 'last'), ...
        find(controls == 2, 1, 'last'), numel(states)];
    assert(all(snapshotIndices > 0), ...
        'The comparison requires injection, storage, and production periods.');
    snapshotLabels = ["End injection", "End storage", "End simulation"];

    aqueousH2 = zeros(model.G.cells.num, numel(snapshotIndices));
    aqueousCO2 = zeros(model.G.cells.num, numel(snapshotIndices));
    for i = 1:numel(snapshotIndices)
        state = states{snapshotIndices(i)};
        aqueousH2(:, i) = componentColumn(state.x, idxH2);
        aqueousCO2(:, i) = componentColumn(state.x, idxCO2);
    end

    nReactions = model.biochemFluid.nbioreact;
    cumulative = zeros(numel(states), nReactions);
    finalSpatial = zeros(model.G.cells.num, nReactions);
    for reaction = 1:nReactions
        [~, reactionCumulative] = computeH2Consumption( ...
            states, schedule, model, reaction);
        cumulative(:, reaction) = sum(reactionCumulative, 1).';
        finalSpatial(:, reaction) = reactionCumulative(:, end);
    end
    totalConsumed = sum(cumulative, 2);
    injected = prescribedInjectedH2(schedule, model);
    x = model.G.cells.centroids(:, 1);
    x = (x - min(x))./max(max(x) - min(x), eps);

    metrics = struct( ...
        'timeDays', cumsum(schedule.step.val(:))./day, ...
        'snapshotIndices', snapshotIndices, ...
        'snapshotLabels', snapshotLabels, ...
        'aqueousH2', aqueousH2, ...
        'aqueousCO2', aqueousCO2, ...
        'lossPercent', 100.*totalConsumed./injected, ...
        'finalSpatialConsumptionMoles', sum(finalSpatial, 2), ...
        'xDimensionless', x, ...
        'injectedH2Moles', injected, ...
        'consumedH2Moles', totalConsumed(end));
end

function scenario = runUGFACTReference(ugfactRoot, useSoreideWhitsonEOS)
    referenceFile = fullfile(ugfactRoot, 'examples', 'H2Storage1D.m');
    assert(isfile(referenceFile), ...
        'UGFACT reference driver was not found at %s.', referenceFile);
    mrstPath('register', 'UGFACT2', ugfactRoot);

    source = fileread(referenceFile);
    assert(contains(source, 'rateCase = ''highrate'''), ...
        'UGFACT H2Storage1D must explicitly select the high-rate case.');
    eosExpression = 'useSoreideWhitsonEOS\s*=\s*(true|false)\s*;';
    assert(~isempty(regexp(source, eosExpression, 'once')), ...
        ['The UGFACT driver EOS selector changed. Update the reference ', ...
         'override before running this comparison.']);
    source = strrep(source, 'clear;clc;close all', '');
    eosValue = char(string(useSoreideWhitsonEOS));
    source = regexprep(source, eosExpression, ...
        ['useSoreideWhitsonEOS = ', eosValue, ';'], 'once');

    if useSoreideWhitsonEOS
        eosName = "Soreide-Whitson EOS";
    else
        eosName = "original UGFACT EOS";
    end
    fprintf('\nBackend 4/4: UGFACT H2Storage1D reference (high rate, %s)\n', ...
        eosName);
    timer = tic();
    evalin('base', source);
    rateCase = evalin('base', 'rateCase');
    model = evalin('base', 'model');
    schedule = evalin('base', 'schedule');
    states = evalin('base', 'states');
    state0 = evalin('base', 'state0');
    wellSol = evalin('base', 'wellSol');
    assert(strcmp(rateCase, 'highrate'), ...
        'UGFACT reference did not run the high-rate kinetic case.');
    assert(model.Kinetic.mu_MET == 4.1 && ...
        model.Kinetic.mu_ACE == 1.9 && model.Kinetic.mu_SRB == 5.5, ...
        'UGFACT reference high-rate kinetic constants are incorrect.');
    actualEOSChoice = evalin('base', 'useSoreideWhitsonEOS');
    assert(logical(actualEOSChoice) == useSoreideWhitsonEOS, ...
        'UGFACT reference did not use the requested EOS.');

    metrics = collectUGFACTMetrics(states, schedule, model);
    scenario = initializeScenarios(1);
    scenario.name = "UGFACT reference (high rate, " + eosName + ")";
    scenario.model = model;
    scenario.schedule = schedule;
    scenario.states = states;
    scenario.ws = wellSol;
    scenario.state0 = state0;
    scenario.timeDays = metrics.timeDays;
    scenario.snapshotIndices = metrics.snapshotIndices;
    scenario.snapshotLabels = metrics.snapshotLabels;
    scenario.aqueousH2 = metrics.aqueousH2;
    scenario.aqueousCO2 = metrics.aqueousCO2;
    scenario.lossPercent = metrics.lossPercent;
    scenario.finalSpatialConsumptionMoles = ...
        metrics.finalSpatialConsumptionMoles;
    scenario.xDimensionless = metrics.xDimensionless;
    scenario.injectedH2Moles = metrics.injectedH2Moles;
    scenario.consumedH2Moles = metrics.consumedH2Moles;
    scenario.runtimeSeconds = toc(timer);
end

function metrics = collectUGFACTMetrics(states, schedule, model)
    assert(numel(states) == numel(schedule.step.val), ...
        'UGFACT states and schedule must contain the same number of steps.');
    componentNames = model.EOSModel.getComponentNames();
    idxH2 = findComponentIndex(componentNames, {'H2', 'Hydrogen'});
    idxCO2 = findComponentIndex(componentNames, {'CO2', 'CarbonDioxide'});

    controls = schedule.step.control(:);
    snapshotIndices = [find(controls == 1, 1, 'last'), ...
        find(controls == 2, 1, 'last'), numel(states)];
    snapshotLabels = ["End injection", "End storage", "End simulation"];
    aqueousH2 = zeros(model.G.cells.num, 3);
    aqueousCO2 = zeros(model.G.cells.num, 3);
    for i = 1:3
        state = states{snapshotIndices(i)};
        aqueousH2(:, i) = componentColumn(state.x, idxH2);
        aqueousCO2(:, i) = componentColumn(state.x, idxCO2);
    end

    cumulative = zeros(numel(states), 1);
    finalSpatial = zeros(model.G.cells.num, 1);
    for step = 1:numel(states)
        state = states{step};
        solution = state.Solution;
        assert(all(isfield(solution, ...
            {'MET_Rate', 'ACE_Rate', 'SRB_Rate', 'Water'})), ...
            'UGFACT state %d is missing reaction-rate diagnostics.', step);
        rateMolesPerDay = (solution.MET_Rate + solution.ACE_Rate + ...
            solution.SRB_Rate).*solution.Water./1000;
        increment = rateMolesPerDay.*schedule.step.val(step)./day;
        finalSpatial = finalSpatial + increment;
        cumulative(step) = sum(finalSpatial);
    end

    injected = prescribedInjectedH2(schedule, model);
    x = model.G.cells.centroids(:, 1);
    x = (x - min(x))./max(max(x) - min(x), eps);
    metrics = struct( ...
        'timeDays', cumsum(schedule.step.val(:))./day, ...
        'snapshotIndices', snapshotIndices, ...
        'snapshotLabels', snapshotLabels, ...
        'aqueousH2', aqueousH2, ...
        'aqueousCO2', aqueousCO2, ...
        'lossPercent', 100.*cumulative./injected, ...
        'finalSpatialConsumptionMoles', finalSpatial, ...
        'xDimensionless', x, ...
        'injectedH2Moles', injected, ...
        'consumedH2Moles', cumulative(end));
end

function column = componentColumn(composition, componentIndex)
    if iscell(composition)
        column = value(composition{componentIndex});
    else
        column = value(composition(:, componentIndex));
    end
    column = column(:);
end

function injected = prescribedInjectedH2(schedule, model)
    names = model.EOSModel.CompositionalMixture.names;
    idxH2 = findComponentIndex(names, {'H2', 'Hydrogen'});
    gasIndex = model.getVaporIndex();
    injected = 0;
    for step = 1:numel(schedule.step.val)
        wells = schedule.control(schedule.step.control(step)).W;
        for w = 1:numel(wells)
            if wells(w).sign <= 0 || ...
                    (isfield(wells, 'status') && ~wells(w).status)
                continue;
            end

            if strcmpi(wells(w).type, 'grat')
                gasRate = wells(w).val;
            elseif strcmpi(wells(w).type, 'rate')
                gasRate = wells(w).val*wells(w).compi(gasIndex);
            else
                continue;
            end
            molarRate = gasRate*model.FacilityModel.pressure/ ...
                (8.314462618*model.FacilityModel.T);
            injected = injected + molarRate*wells(w).components(idxH2)* ...
                schedule.step.val(step);
        end
    end
    assert(injected > 0, 'No prescribed H2 injection was found.');
end

function index = findComponentIndex(names, aliases)
    for i = 1:numel(aliases)
        candidate = find(strcmpi(names, aliases{i}), 1);
        if ~isempty(candidate)
            index = candidate;
            return;
        end
    end
    error('Component "%s" was not found.', strjoin(aliases, '" or "'));
end

function printSummary(scenarios)
    names = [scenarios.name].';
    injected = [scenarios.injectedH2Moles].';
    consumed = [scenarios.consumedH2Moles].';
    loss = 100.*consumed./injected;
    runtime = [scenarios.runtimeSeconds].';
    summary = table(names, injected, consumed, loss, runtime, ...
        'VariableNames', {'Backend', 'InjectedH2_mol', 'ConsumedH2_mol', ...
        'FinalH2Loss_percent', 'Runtime_seconds'});
    disp(summary);
end

function createComparisonFigures(scenarios)
    colors = [ ...
        0.15, 0.15, 0.15; ...
        0.00, 0.45, 0.74; ...
        0.85, 0.33, 0.10; ...
        0.47, 0.67, 0.19];
    labels = scenarios(1).snapshotLabels;

    fig = paperFigure([24, 15], 'Aqueous H2 and CO2 profiles');
    layout = tiledlayout(fig, 2, 3, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    for component = 1:2
        for snapshot = 1:3
            ax = nexttile(layout);
            hold(ax, 'on');
            for backend = 1:numel(scenarios)
                if component == 1
                    profile = scenarios(backend).aqueousH2(:, snapshot);
                else
                    profile = scenarios(backend).aqueousCO2(:, snapshot);
                end
                plot(ax, scenarios(backend).xDimensionless, profile, ...
                    'LineWidth', 1.5, 'Color', colors(backend, :), ...
                    'DisplayName', scenarios(backend).name);
            end
            title(ax, labels(snapshot));
            if component == 1
                ylabel(ax, 'Aqueous H_2 mole fraction');
            else
                ylabel(ax, 'Aqueous CO_2 mole fraction');
            end
            xlabel(ax, 'Dimensionless distance');
            styleAxes(ax);
            if component == 1 && snapshot == 3
                legend(ax, 'Location', 'best', 'Box', 'off');
            end
        end
    end

    fig = paperFigure([18, 11], 'H2 loss over time');
    ax = axes(fig);
    hold(ax, 'on');
    for backend = 1:numel(scenarios)
        plot(ax, scenarios(backend).timeDays, ...
            scenarios(backend).lossPercent, ...
            'LineWidth', 1.6, 'Color', colors(backend, :), ...
            'DisplayName', scenarios(backend).name);
    end
    xline(ax, 50, 'k:', 'End injection');
    xline(ax, 200, 'k:', 'End storage');
    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Consumed injected H_2 (%)');
    legend(ax, 'Location', 'best', 'Box', 'off');
    styleAxes(ax);

    fig = paperFigure([18, 11], 'Spatial H2 consumption');
    ax = axes(fig);
    hold(ax, 'on');
    for backend = 1:numel(scenarios)
        plot(ax, scenarios(backend).xDimensionless, ...
            scenarios(backend).finalSpatialConsumptionMoles, ...
            'LineWidth', 1.6, 'Color', colors(backend, :), ...
            'DisplayName', scenarios(backend).name);
    end
    xlabel(ax, 'Dimensionless distance from injector');
    ylabel(ax, 'Cumulative H_2 consumed (mol cell^{-1})');
    legend(ax, 'Location', 'best', 'Box', 'off');
    styleAxes(ax);
end

function fig = paperFigure(sizeCm, name)
    fig = figure('Name', name, 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2, 2, sizeCm]);
end

function styleAxes(ax)
    grid(ax, 'on');
    box(ax, 'on');
    ax.FontName = 'Times New Roman';
    ax.FontSize = 9;
    ax.LineWidth = 0.8;
    ax.TickDir = 'out';
    ax.Layer = 'top';
end

function databaseFile = resolvePhreeqcDatabaseFile(databaseFile)
    assert(ischar(databaseFile) || ...
        (isstring(databaseFile) && isscalar(databaseFile)), ...
        'phreeqcDatabaseFile must be a character vector or scalar string.');
    databaseFile = char(databaseFile);
    if isempty(strtrim(databaseFile))
        databaseFile = getenv('PHREEQC_DATABASE_FILE');
    end
    if isempty(strtrim(databaseFile))
        databaseFile = which('PHREEQC_Modified.DAT');
    end
    if isempty(strtrim(databaseFile)) || ~isfile(databaseFile)
        error('runThreeBackendComparison:MissingPhreeqcDatabase', ...
            ['PHREEQC_Modified.DAT was not found. Pass its absolute path ', ...
            'using ''phreeqcDatabaseFile'', databaseFile.']);
    end
    [~, name, extension] = fileparts(databaseFile);
    assert(strcmpi([name, extension], 'PHREEQC_Modified.DAT'), ...
        'The database must be PHREEQC_Modified.DAT.');
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
%}
