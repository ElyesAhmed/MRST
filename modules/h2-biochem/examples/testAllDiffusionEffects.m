function scenarios = testAllDiffusionEffects()
% Test all diffusion/dispersion/chemotaxis/bio-clogging effects in h2-biochem module.
% Runs setupH2StorageExample with all combinations of flags and compares results.
%
% DESCRIPTION:
%   This script tests the h2-biochem module by running the setupH2StorageExample case
%   with different combinations of transport physics flags:
%   - bacteriamodel: master switch for bacterial growth/decay
%   - bactDiffusion: microbial (Fickian) diffusion
%   - chemotaxisEffect: bacterial chemotaxis
%   - molecularDiffusion: molecular diffusion of components
%   - molecularDispersion: mechanical dispersion
%   - bioClogging: porosity/permeability reduction from bacterial activity
%
%   Results from each run are collected into a scenarios struct and plotted
%   using existing utilities (plotBacterialEffects, plotComponentProfiles, plotH2Loss).

    mrstModule add ad-props compositional h2-biochem

    fprintf('\n===== Testing all diffusion effects in h2-biochem =====\n\n');

    % Define flag combinations to test
    % Start with baseline (all off), then enable each effect individually,
    % then enable combinations
    flagSets = {
        struct('bacteriamodel', false, 'bactDiffusion', false, 'chemotaxisEffect', false, ...
               'molecularDiffusion', false, 'molecularDispersion', false, 'bioClogging', false,...
               'phreeqcTimestepCoupling', false, 'phreeqcConservativeCarbonTransfer', false, 'phreeqcInitialization', false, 'carbonateBuffer', false), ...
        struct('bacteriamodel', true, 'bactDiffusion', false, 'chemotaxisEffect', false, ...
               'molecularDiffusion', false, 'molecularDispersion', false, 'bioClogging', false, ...
               'phreeqcTimestepCoupling', true, 'phreeqcConservativeCarbonTransfer', true, 'phreeqcInitialization', true, 'carbonateBuffer', true), ...
        struct('bacteriamodel', true, 'bactDiffusion', true, 'chemotaxisEffect', false, ...
               'molecularDiffusion', false, 'molecularDispersion', false, 'bioClogging', false), ...
        struct('bacteriamodel', true, 'bactDiffusion', false, 'chemotaxisEffect', true, ...
               'molecularDiffusion', false, 'molecularDispersion', false, 'bioClogging', false), ...
        struct('bacteriamodel', true, 'bactDiffusion', true, 'chemotaxisEffect', true, ...
               'molecularDiffusion', false, 'molecularDispersion', false, 'bioClogging', false), ...
        struct('bacteriamodel', false, 'bactDiffusion', false, 'chemotaxisEffect', false, ...
               'molecularDiffusion', true, 'molecularDispersion', false, 'bioClogging', false), ...
        struct('bacteriamodel', false, 'bactDiffusion', false, 'chemotaxisEffect', false, ...
               'molecularDiffusion', false, 'molecularDispersion', true, 'bioClogging', false), ...
        struct('bacteriamodel', false, 'bactDiffusion', false, 'chemotaxisEffect', false, ...
               'molecularDiffusion', true, 'molecularDispersion', true, 'bioClogging', false), ...
        struct('bacteriamodel', true, 'bactDiffusion', true, 'chemotaxisEffect', false, ...
               'molecularDiffusion', true, 'molecularDispersion', true, 'bioClogging', false), ...
        struct('bacteriamodel', true, 'bactDiffusion', true, 'chemotaxisEffect', true, ...
               'molecularDiffusion', true, 'molecularDispersion', true, 'bioClogging', true), ...
    };

    % Scenario names for plotting
    sceneNames = {
        'Baseline (all off)', ...
        'Bacteria only', ...
        'Bacteria + Diffusion', ...
        'Bacteria + Chemotaxis', ...
        'Bacteria + Diffusion + Chemotaxis', ...
        'Molecular Diffusion', ...
        'Mechanical Dispersion', ...
        'Dispersion + Diffusion', ...
        'Full Physics (no clogging)', ...
        'Full Physics (with clogging)', ...
    };

    % Color map for scenarios
    colors = lines(length(flagSets));

    scenarios = {};
    figNum = 0;

    % Run each test case
    for i = 2:length(flagSets)
        flags = flagSets{i};
        fprintf('Test %d/%d: %s\n', i, length(flagSets), sceneNames{i});

        %
            % Set up model with this flag combination
            [biochemFluid, model, schedule, state0] = setupH2StorageExampleWithSRB_benchmark( ...
                'phreeqcInitialization', flags.phreeqcInitialization, ...
                'phreeqcTimestepCoupling', flags.phreeqcTimestepCoupling, ...
                'phreeqcConservativeCarbonTransfer', flags.phreeqcConservativeCarbonTransfer, ...
                'injectionCO2', 0, ...
                'bacteriamodel', flags.bacteriamodel, ...
                'bactDiffusion', flags.bactDiffusion, ...
                'chemotaxisEffect', flags.chemotaxisEffect, ...
                'molecularDiffusion', flags.molecularDiffusion, ...
                'molecularDispersion', flags.molecularDispersion, ...
                'bioClogging', flags.bioClogging);

            % Run simulation
            fprintf('  Running simulation...');
            [ws, states] = simulateScheduleAD(state0, model, schedule);
            fprintf(' OK\n');

            % Store results in scenarios struct
            scenario = struct();
            scenario.name = sceneNames{i};
            scenario.state0 = state0;
            scenario.schedule = schedule;
            scenario.states = states;
            scenario.ws = ws;
            scenario.model = model;
            scenario.color = colors(i, :);
            scenario.line = '-';

            scenarios{end+1} = scenario;

    end

    if isempty(scenarios)
        error('No scenarios completed successfully');
    end

    fprintf('\nCompleted %d/%d test cases successfully\n\n', length(scenarios), length(flagSets));

    % Plot results
    fprintf('Generating plots...\n');

    % Extract simulation time in years
    dt = schedule.step.val;
    timeYears = cumsum(dt) / year;

    % Plot 1: Bacterial effects (if any bacterial simulations)
    hasBacteria = any(cellfun(@(s) s.model.bacteriamodel, scenarios));
    if hasBacteria
        figNum = figNum + 1;
        figure(figNum);
        try
            plotBacterialEffects(scenarios, timeYears);
            title('Bacterial Concentration and Porosity/Permeability Changes');
        catch
            fprintf('  Could not plot bacterial effects\n');
        end
    end

    % Plot 2: Component profiles
    figNum = figNum + 1;
    figure(figNum);
    try
        compNames = scenarios{1}.model.compFluid.names;
        plotComponentProfiles(scenarios, compNames, timeYears);
        title('Component Composition Changes');
    catch
        fprintf('  Could not plot component profiles\n');
    end

    % Plot 3: H2 loss
    figNum = figNum + 1;
    figure(figNum);
    try
        plotH2Loss(scenarios{1}.model, schedule, states, ws);
        title('H2 Mass Loss Analysis');
    catch
        fprintf('  Could not plot H2 loss\n');
    end

    fprintf('Plots complete.\n\n');
    fprintf('===== Testing finished =====\n');

end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
