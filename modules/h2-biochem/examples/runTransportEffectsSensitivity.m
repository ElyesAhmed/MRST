function results = runTransportEffectsSensitivity(varargin)
% Run a transport and bio-clogging sensitivity study and display its plots.
%
% The baseline is the first "Original h2-biochem" case from
% testAllDiffusionEffects: MRST owns all reactions, the reduced fixed-pH
% CO2/HCO3 buffer is enabled, and PHREEQC timestep coupling is disabled.
% Results remain in memory and figures are left open; no files are written.

    mrstModule add ad-props compositional deckformat h2-biochem

    opt = struct( ...
        'meaningfulLossThreshold', 0.1, ...
        'scheduleMode', 'complete', ...
        'rate', 'highrate');
    opt = merge_options(opt, varargin{:});
    validateOptions(opt);

    cases = transportCases();
    results = initializeResults(numel(cases));

    commonOptions = { ...
        'rate', opt.rate, ...
        'paperBiomassKinetics', false, ...
        'scheduleMode', opt.scheduleMode, ...
        'injectionCO2', 0, ...
        'carbonateBuffer', true, ...
        'equilibrateInitialCO2', true, ...
        'phreeqcTimestepCoupling', false, ...
        'bacteriamodel', true};

    for caseNo = 1:numel(cases)
        fprintf('\nCase %d/%d: %s\n', caseNo, numel(cases), ...
            cases(caseNo).name);
        setupOptions = [commonOptions, cases(caseNo).options];
        [~, model, schedule, state0] = ...
            setupH2StorageExampleWithSRB_benchmark(setupOptions{:});
        assert(~model.phreeqcTimestepCoupling, ...
            'Transport sensitivity cases must not call PHREEQC.');
        assert(~model.isSequentialCompositionalPhreeqcBackend() && ...
            ~model.isSequentialH2BiochemPhreeqcBackend(), ...
            'Transport sensitivity cases must use MRST-owned reactions.');
        solver = NonLinearSolver();
        solver.maxTimestepCuts = 12;
        timer = tic();
        [~, states, report] = simulateScheduleAD( ...
            state0, model, schedule, 'nonlinearSolver', solver);
        assert(~report.Failure, ...
            'Simulation report indicates failure for case "%s".', ...
            cases(caseNo).name);
        metrics = computeMetrics(states, schedule, model);
        results(caseNo).completed = true;
        results(caseNo).runtimeSeconds = toc(timer);
        results(caseNo).timeDays = metrics.timeDays;
        results(caseNo).lossPercent = metrics.lossPercent;
        results(caseNo).finalLossPercent = metrics.finalLossPercent;
        results(caseNo).consumedH2Moles = metrics.consumedH2Moles;
        results(caseNo).injectedH2Moles = metrics.injectedH2Moles;
        results(caseNo).reactionTotalsMoles = metrics.reactionTotalsMoles;
        results(caseNo).finalSpatialConsumptionMoles = ...
            metrics.finalSpatialConsumptionMoles;
        results(caseNo).xDimensionless = metrics.xDimensionless;
        results(caseNo).errorMessage = "";
    end

    summary = buildSummary(results, cases, opt);
    createPaperFigures(results, cases, summary, opt);
    disp(summary);
end

function validateOptions(opt)
    validateattributes(opt.meaningfulLossThreshold, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'}, ...
        mfilename, 'meaningfulLossThreshold');
end

function cases = transportCases()
    cases = [ ...
        makeCase("Baseline", false, false, false, false, false), ...
        makeCase("Microbial diffusion", true, false, false, false, false), ...
        makeCase("Chemotaxis", false, true, false, false, false), ...
        makeCase("Microbial diffusion + chemotaxis", ...
            true, true, false, false, false), ...
        makeCase("Molecular diffusion", false, false, true, false, false), ...
        makeCase("Mechanical dispersion", false, false, false, true, false), ...
        makeCase("Molecular diffusion + dispersion", ...
            false, false, true, true, false), ...
        makeCase("Bio-clogging", false, false, false, false, true), ...
        makeCase("All transport effects", true, true, true, true, false), ...
        makeCase("All transport effects + bio-clogging", ...
            true, true, true, true, true)];
end

function item = makeCase(name, bactDiffusion, chemotaxis, ...
        molecularDiffusion, mechanicalDispersion, bioClogging)
    item = struct( ...
        'name', char(name), ...
        'options', {{ ...
        'bactDiffusion', bactDiffusion, ...
        'chemotaxisEffect', chemotaxis, ...
        'molecularDiffusion', molecularDiffusion, ...
        'molecularDispersion', mechanicalDispersion, ...
        'bioClogging', bioClogging}});
end

function results = initializeResults(nCases)
    template = struct( ...
        'completed', false, ...
        'runtimeSeconds', nan, ...
        'timeDays', [], ...
        'lossPercent', [], ...
        'finalLossPercent', nan, ...
        'consumedH2Moles', nan, ...
        'injectedH2Moles', nan, ...
        'reactionTotalsMoles', [], ...
        'finalSpatialConsumptionMoles', [], ...
        'xDimensionless', [], ...
        'errorMessage', "");
    results = repmat(template, nCases, 1);
end

function metrics = computeMetrics(states, schedule, model)
    nReactions = model.biochemFluid.nbioreact;
    cumulative = zeros(numel(states), nReactions);
    finalSpatial = zeros(model.G.cells.num, nReactions);
    for reaction = 1:nReactions
        [~, reactionCumulative] = computeH2Consumption( ...
            states, schedule, model, reaction);
        cumulative(:, reaction) = sum(reactionCumulative, 1).';
        finalSpatial(:, reaction) = reactionCumulative(:, end);
    end
    total = sum(cumulative, 2);
    injected = prescribedInjectedH2(schedule, model);
    x = model.G.cells.centroids(:, 1);
    x = (x - min(x))./max(max(x) - min(x), eps);
    metrics = struct( ...
        'timeDays', cumsum(schedule.step.val(:))./day, ...
        'lossPercent', 100.*total./injected, ...
        'finalLossPercent', 100.*total(end)./injected, ...
        'consumedH2Moles', total(end), ...
        'injectedH2Moles', injected, ...
        'reactionTotalsMoles', cumulative(end, :), ...
        'finalSpatialConsumptionMoles', sum(finalSpatial, 2), ...
        'xDimensionless', x);
end

function injected = prescribedInjectedH2(schedule, model)
    names = model.EOSModel.CompositionalMixture.names;
    idxH2 = find(strcmpi(names, 'H2'), 1);
    assert(~isempty(idxH2), 'The compositional mixture does not contain H2.');
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

function summary = buildSummary(results, cases, opt)
    finalLoss = [results.finalLossPercent].';
    delta = finalLoss - finalLoss(1);
    relative = 100.*delta./max(abs(finalLoss(1)), eps);
    direction = strings(numel(cases), 1);
    meaningful = strings(numel(cases), 1);
    for i = 1:numel(cases)
        if abs(delta(i)) < opt.meaningfulLossThreshold
            direction(i) = "No material change";
            meaningful(i) = "No";
        elseif delta(i) > 0
            direction(i) = "Increased loss";
            meaningful(i) = "Yes";
        else
            direction(i) = "Decreased loss";
            meaningful(i) = "Yes";
        end
    end
    summary = table(string({cases.name}).', finalLoss, delta, relative, ...
        [results.consumedH2Moles].', [results.runtimeSeconds].', ...
        direction, meaningful, ...
        'VariableNames', {'Case', 'FinalH2Loss_percent', ...
        'DeltaFromBaseline_percentagePoints', 'RelativeChange_percent', ...
        'ConsumedH2_mol', 'Runtime_seconds', 'Direction', 'Meaningful'});
end

function createPaperFigures(results, cases, summary, opt)
    colors = lines(numel(results));
    fig = paperFigure([18, 11]);
    ax = axes(fig);
    hold(ax, 'on');
    for i = 1:numel(results)
        plot(ax, results(i).timeDays, results(i).lossPercent, ...
            'LineWidth', 1.35, 'Color', colors(i, :), ...
            'DisplayName', cases(i).name);
    end
    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Consumed injected H_2 (%)');
    legend(ax, 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 7);
    styleAxes(ax);

    fig = paperFigure([18, 12]);
    ax = axes(fig);
    values = summary.DeltaFromBaseline_percentagePoints;
    labels = categorical(summary.Case, summary.Case);
    bars = barh(ax, labels, values, 'FaceColor', 'flat');
    bars.CData = repmat([0.16, 0.48, 0.70], numel(values), 1);
    bars.CData(values < 0, :) = repmat([0.78, 0.30, 0.22], ...
        nnz(values < 0), 1);
    xline(ax, 0, 'k-', 'LineWidth', 0.8);
    xline(ax, opt.meaningfulLossThreshold, 'k:', 'LineWidth', 1);
    xline(ax, -opt.meaningfulLossThreshold, 'k:', 'LineWidth', 1);
    xlabel(ax, '\Delta final H_2 loss (percentage points)');
    styleAxes(ax);

    selected = [1, 2, 3, 5, 6, 9, 10];
    fig = paperFigure([18, 11]);
    ax = axes(fig);
    hold(ax, 'on');
    for i = selected
        plot(ax, results(i).xDimensionless, ...
            results(i).finalSpatialConsumptionMoles, ...
            'LineWidth', 1.35, 'Color', colors(i, :), ...
            'DisplayName', cases(i).name);
    end
    xlabel(ax, 'Dimensionless distance from injector');
    ylabel(ax, 'Cumulative H_2 consumption (mol cell^{-1})');
    legend(ax, 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 7);
    styleAxes(ax);
end

function fig = paperFigure(sizeCm)
    fig = figure('Color', 'w', 'Units', 'centimeters', ...
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

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
%}
