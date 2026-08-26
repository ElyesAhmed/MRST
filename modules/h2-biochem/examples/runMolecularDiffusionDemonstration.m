function results = runMolecularDiffusionDemonstration()
% Compare H2 loss with and without molecular diffusion in a lab-scale column.
%
% The shorter column and long shut-in make molecular transport resolvable
% without changing the physical diffusion coefficients. No files are written.

    mrstModule add ad-props compositional deckformat h2-biochem

    caseNames = ["No molecular diffusion"; "Molecular diffusion"];
    diffusionEnabled = [false; true];
    nCases = numel(caseNames);
    results = repmat(struct( ...
        'name', "", ...
        'runtimeSeconds', nan, ...
        'timeDays', [], ...
        'lossPercent', [], ...
        'consumedH2Moles', nan, ...
        'injectedH2Moles', nan, ...
        'finalSpatialConsumptionMoles', [], ...
        'distanceMetres', []), nCases, 1);

    commonOptions = { ...
        'rate', 'highrate', ...
        'paperBiomassKinetics', false, ...
        'scheduleMode', 'complete', ...
        'gridCells', 50, ...
        'domainLength', 10, ...
        'injectionCO2', 0, ...
        'carbonateBuffer', true, ...
        'equilibrateInitialCO2', true, ...
        'phreeqcTimestepCoupling', false, ...
        'bacteriamodel', true, ...
        'bactDiffusion', false, ...
        'chemotaxisEffect', false, ...
        'molecularDispersion', false, ...
        'bioClogging', false};

    for caseNo = 1:nCases
        fprintf('\nMolecular-diffusion case %d/%d: %s\n', ...
            caseNo, nCases, caseNames(caseNo));
        [~, model, schedule, state0] = ...
            setupH2StorageExampleWithSRB_benchmark(commonOptions{:}, ...
            'molecularDiffusion', diffusionEnabled(caseNo));
        assert(~model.phreeqcTimestepCoupling, ...
            'The molecular-diffusion demonstration must not call PHREEQC.');
        schedule = extendStoragePeriod(schedule, 500*day, 5*day);

        solver = NonLinearSolver();
        solver.maxTimestepCuts = 12;
        timer = tic();
        [~, states, report] = simulateScheduleAD( ...
            state0, model, schedule, 'nonlinearSolver', solver);
        assert(~report.Failure, ...
            'Simulation report indicates failure for case "%s".', ...
            caseNames(caseNo));

        metrics = computeMetrics(states, schedule, model);
        results(caseNo).name = caseNames(caseNo);
        results(caseNo).runtimeSeconds = toc(timer);
        results(caseNo).timeDays = metrics.timeDays;
        results(caseNo).lossPercent = metrics.lossPercent;
        results(caseNo).consumedH2Moles = metrics.consumedH2Moles;
        results(caseNo).injectedH2Moles = metrics.injectedH2Moles;
        results(caseNo).finalSpatialConsumptionMoles = ...
            metrics.finalSpatialConsumptionMoles;
        results(caseNo).distanceMetres = model.G.cells.centroids(:, 1);
    end

    finalLoss = arrayfun(@(r) r.lossPercent(end), results);
    deltaLoss = finalLoss - finalLoss(1);
    summary = table(caseNames, finalLoss, deltaLoss, ...
        [results.consumedH2Moles].', [results.runtimeSeconds].', ...
        'VariableNames', {'Case', 'FinalH2Loss_percent', ...
        'DeltaFromBaseline_percentagePoints', 'ConsumedH2_mol', ...
        'Runtime_seconds'});
    disp(summary);
    createFigures(results);
end

function schedule = extendStoragePeriod(schedule, duration, maxStep)
    injection = schedule.step.control == 1;
    production = schedule.step.control == 3;
    storageSteps = repmat(maxStep, floor(duration/maxStep), 1);
    remainder = duration - sum(storageSteps);
    if remainder > 0
        storageSteps(end + 1, 1) = remainder;
    end
    schedule.step.val = [schedule.step.val(injection); ...
        storageSteps; schedule.step.val(production)];
    schedule.step.control = [ones(nnz(injection), 1); ...
        2*ones(numel(storageSteps), 1); ...
        3*ones(nnz(production), 1)];
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
    metrics = struct( ...
        'timeDays', cumsum(schedule.step.val(:))./day, ...
        'lossPercent', 100.*total./injected, ...
        'consumedH2Moles', total(end), ...
        'injectedH2Moles', injected, ...
        'finalSpatialConsumptionMoles', sum(finalSpatial, 2));
end

function injected = prescribedInjectedH2(schedule, model)
    names = model.EOSModel.CompositionalMixture.names;
    idxH2 = find(strcmpi(names, 'H2'), 1);
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

function createFigures(results)
    colors = [0.25, 0.25, 0.25; 0.00, 0.45, 0.74];

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2, 2, 18, 15]);
    layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax = nexttile(layout);
    hold(ax, 'on');
    for i = 1:numel(results)
        plot(ax, results(i).timeDays, results(i).lossPercent, ...
            'LineWidth', 1.6, 'Color', colors(i, :), ...
            'DisplayName', results(i).name);
    end
    xline(ax, 50, 'k:', 'Injection ends', ...
        'LabelVerticalAlignment', 'bottom');
    xline(ax, 550, 'k:', 'Production begins', ...
        'LabelVerticalAlignment', 'bottom');
    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Consumed injected H_2 (%)');
    legend(ax, 'Location', 'best', 'Box', 'off');
    styleAxes(ax);

    ax = nexttile(layout);
    delta = results(2).lossPercent - results(1).lossPercent;
    plot(ax, results(1).timeDays, delta, 'LineWidth', 1.6, ...
        'Color', colors(2, :));
    yline(ax, 0, 'k-');
    xline(ax, 50, 'k:');
    xline(ax, 550, 'k:');
    xlabel(ax, 'Time (days)');
    ylabel(ax, '\Delta H_2 loss (percentage points)');
    styleAxes(ax);

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2, 2, 18, 11]);
    ax = axes(fig);
    hold(ax, 'on');
    for i = 1:numel(results)
        plot(ax, results(i).distanceMetres, ...
            results(i).finalSpatialConsumptionMoles, ...
            'LineWidth', 1.6, 'Color', colors(i, :), ...
            'DisplayName', results(i).name);
    end
    xlabel(ax, 'Distance from injector (m)');
    ylabel(ax, 'Cumulative H_2 consumption (mol cell^{-1})');
    legend(ax, 'Location', 'best', 'Box', 'off');
    styleAxes(ax);
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
