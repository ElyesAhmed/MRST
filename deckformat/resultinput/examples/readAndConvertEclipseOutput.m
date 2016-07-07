%% Read ECLIPSE output and convert to MRST grid and solution structures
% This example reads ECLIPSE unformatted output files from a simulation based 
% on the SPE9 benchmark. The MRST grid structure (G) is generated using output
% from *.INIT and *.EGRID. The solution structure (states and states{i}.wellSol)
% are generated using output from *.UNRST.

mrstModule add mrst-gui ad-core
prefix =  '/Users/steink/data/opm-data/spe9/eclipse-simulation/SPE9_CP';

%% Read INIT/EGRID-files and generate MRST-grid
init =  readEclipseOutputFileUnFmt([prefix, '.INIT']);
egrid = readEclipseOutputFileUnFmt([prefix, '.EGRID']);
[G, rock] = eclOut2mrst(init, egrid);
G = computeGeometry(G);

%% Restart-to-states converter reads selected parts of UNRST-file.
% Each state (states{i}) contains a 'wellSol' structure containing
% sufficient welldata for plotting (e.g., perforation cells). 
states = convertRestartToStates(prefix, G);

% Fetch wellsols and restart-step times
wellSols = cellfun(@(x)x.wellSol, states, 'UniformOutput', false);
time     = cellfun(@(x)x.time, states);

%% Plot model with wells
figure, plotToolbar(G, states);
plotWellData(G, wellSols{1}, 'labelFontSize', 10, 'color', [.3 .3 .3]);
axis off vis3d tight, view([1 -.5 1]), camproj perspective

%% Plot well solutions for restart time-steps
plotWellSols(wellSols, time);

