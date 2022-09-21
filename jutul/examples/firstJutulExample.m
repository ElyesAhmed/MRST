%% First introductory example to Jutul as a MRST accelerator
% JutulDarcy is a reservoir simulator written in Julia by SINTEF Digital
% based on on the Jutul solver framework. We can set up a case in MRST, and
% run it in Jutul for increased computational performance. Jutul supports
% black-oil, immiscible and compositional models with multisegment wells.
%
% To install Jutul, first install the latest version of Julia [1].
% Once installed, run Julia and add the JutulDarcy package . If you are
% interested in using Julia, for other things, we recommend adding it to an
% environment[3]. Otherwise, you can add it to the default environment by
% running the following command in the Julia prompt:
%
% using Pkg; Pkg.add("JutulDarcy")
%
% Once downloaded, you are ready to run this example. Note that the example
% pauses once the simulation is ready to be run in the Julia terminal.
%
% For more details on Jutul, see the JutulDarcy repository on GitHub [3]

% [1] https://julialang.org/downloads/
% [2] https://pkgdocs.julialang.org/v1/environments/
% [3] https://github.com/sintefmath/JutulDarcy.jl

mrstModule add ad-core ad-blackoil spe10 deckformat ad-props test-suite compositional jutul
if ~exist('name', 'var')
    name = 'qfs_wo';
end
%%
switch name
    case 'qfs_wo'
        setup = qfs_wo();
    case 'buckley_leverett_wo'
        setup = buckley_leverett_wo();
    case 'spe10_layer'
        setup = spe10_wo('layers', 1);
    case 'fractures_compositional'
        setup = fractures_compositional();
    otherwise
        % Assume some test-suite case
        setup = eval(name);
end
[state0, model, schedule] = deal(setup.state0, setup.model, setup.schedule);
%% Write case to disk
pth = writeJutulInput(state0, model, schedule, name);
disp('Pausing - run the command in Julia and hit any key to continue')
pause()
%% Once simulated, read back as MRST format
[ws, states] = readJutulOutput(pth);
%% Simulate MRST for comparison purposes
nls = getNonLinearSolver(model);
[ws_m, states_m] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);
%% Compare the results
mrstModule add mrst-gui
G = model.G;
figure;
plotToolbar(G, states_m)
title('MRST')
figure;
plotToolbar(G, states)
title('Jutul')
figure;
plotToolbar(G, applyFunction(@(x, y) compareStructs(x, y), states, states_m))
title('Difference')
%% Plot and compare wells
% Jutul uses multisegment wells by default which gives small differences in
% well curves
plotWellSols({ws, ws_m}, cumsum(schedule.step.val), 'datasetnames', {'Jutul', 'MRST'})