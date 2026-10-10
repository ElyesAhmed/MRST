%% Water injection into a 3D fractured porous media
% Two-phase example with a horizontal producer and injector simulating
% water injection in a 3-dimensional fractured porous media using the HFM
% module. Note that the 3D solvers are not capable of handling intersecting
% fracture planes.

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add ad-core;         % NNC support for coarse grids
mrstModule add msrsb;           % MsRSB solvers
mrstModule add mrst-gui;        % plotting routines
mrstModule add incomp;          % Incompressible fluid models
checkMATLABversionHFM;

%% Grid and fracture(s)
% Construct a Cartesian grid comprising 25-by-25-by-25 cells, where each
% cell has dimension 4-by-4-by-4 m^3. Define a fracture plane using a
% coplanar set of points supplied as rows in the matrix [X Y Z]. In this
% example, we define a fracture plane of size 50-by-50 m^2, extending in
% the y-direction through the centre of the domain. Fracture aperture is
% set to 0.02 meters. Additionally, we compute a triangulation using the
% matrix grid nodes. The triangulation is used to locate cells penetrated
% by the fracture plane in the next section.
celldim = [14, 14, 10];
physdim = [100 100 100];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
GlobTri = globalTriangulation(G);
fracplanes = struct;
fracplanes(1).points = [50 25 25;
                        50 25 75;
                        50 75 75;
                        50 75 25];
fracplanes(1).aperture = 1/50;
checkIfCoplanar(fracplanes)

%% Process fracture(s)
% Using the input fracture, we identify fine-cells in the matrix containing
% these fractures. Next, a global conductivity index (CI) is computed for
% the fracture plane and a grid is defined over it. A 'non-neighboring
% connection (NNC)' status is given to each fracture-matrix connection and
% we compute the transmissibility for each NNC using the global CI.
dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracplanes] = preProcessingFractures(G, fracplanes, ...
                 'GlobTri', GlobTri, ...
                 'fractureCellSize', 0.3);
fprintf('\nProcessing complete...\n');

figure; plotGrid(G,'facealpha',0);
for i = 1:numel(fieldnames(G.FracGrid))
    plotGrid(G.FracGrid.(['Frac',num2str(i)]));
end
view(15,20);

%% Set rock properties in fracture and matrix
% Set the permeability (K) as 1 Darcy in the matrix and 10000 Darcy in the
% fractures. Additionally, set the porosity of the matrix and fractures to
% 20% and 50%, respectively.
dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
rock = getSPE10rock(1:40, 101:149, 4);
G.rock.perm = rock.perm(:,1);%ones(G.cells.num,1)*darcy();
G.rock.poro = rock.poro;
K_frac = 10000; % Darcy
poro_frac = 0.5;
G = makeRockFrac(G, K_frac, 'permtype','homogeneous','porosity', poro_frac);

%% Define fluid properties
% Define a two-phase fluid model without capillarity. The fluid model has
% the values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 1] cP.
% * corey-coefficient: [2, 2] = [2 2].
fluid = initSimpleFluid('mu' , [   1,  1] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

%% Assemble global grid and compute transmissibilities
% In this section, we combine the fracture and matrix grids into one grid.
% The transmissibility for each face in the combined grid and each NNC is
% computed and stored in the vector T.
G = assembleGlobalGrid(G);
G = computeEffectiveTrans(G);
G.faces.tag = zeros(G.faces.num,1);
x = ismember(G.faces.neighbors,G.nnc.cells,'rows');
G.faces.tag(x) = 1; % Tag NNC faces. 
%
T = computeTrans(G, G.rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T = [T;G.nnc.T];

%% Add wells
% We will include two horizontal wells, a rate-controlled injector well
% through the bottom left corner of the grid and a producer, controlled by
% bottom-hole pressure, through the top right section of the grid. Wells
% are described using a Peaceman model, giving an extra set of equations
% that need to be assembled, see simpleWellExample.m for more details on
% the Peaceman well model.
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
cellinj = nx*ny*(nz-1)+1:nx:nx*ny*nz;
cellprod = nx:nx:nx*ny;
W   = addWell([], G.Matrix, G.Matrix.rock, flipud(cellinj), 'Type', 'rate',...
    'Val', 500/day, 'Sign',1, 'Comp_i', [1, 0], 'Name', 'Injector');
W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod, 'Type', 'bhp', ...
    'Val', 100*barsa, 'Sign', -1, 'Comp_i', [0, 1], 'Name', 'Producer');
plotWell(G,W);

%% Set up schedule with 4 control steps
% Following the pattern from original: 4 control steps with varying timesteps
ts = { [1 1 3 5 5 10 10 10 15 15 15 15 15 15 15]'*day, ...
                   repmat(150/10, 10, 1)*day, ...
                   repmat(150/6, 6, 1)*day, ...
                   repmat(150/6, 6, 1)*day};
       
numCnt = numel(ts);
[schedule.control(1:numCnt).W] = deal(W);
schedule.step.control = rldecode((1:4)', cellfun(@numel, ts));
schedule.step.val     = vertcat(ts{:});

gravity off;

% modelRef = GenericBlackOilModel(G.Matrix, G.Matrix.rock, fluid, 'gas', false);
% % Set initial state and run simulation:
% state0 = initResSol(G.Matrix, 200*barsa, [0, 1]);
% % 1.2 Visualize reference model
% figure(1);
% subplot(1,2,1);
% plotCellData(modelRef.G, modelRef.rock.poro, 'EdgeAlpha', 0.5); 
% view(3); axis tight off;
% title('Reference Model Porosity');
% plotWell(modelRef.G, W, 'Color', 'k');
% drawnow;
% 
% % 1.3 Create and perturb schedule for realistic inversion scenario
% time_steps = schedule.step.val;
% W = schedule.control(1).W;
% 
% % Add controlled perturbations to well operations
% rng(0); % For reproducibility
% perturbStep = [ones(1,4), round(.5+(2:(numel(schedule.step.val)-3))/2)]';
% scheduleRef = perturbedSimpleSchedule(time_steps, 'W', W, ...
%     'pressureFac', 0.375, 'rateFac', 0.45, 'perturbStep', perturbStep);
% 
% % 1.4 Configure model for precise simulations
% modelRef.useCNVConvergence = false; 
% modelRef.toleranceMB = 1.0e-9;
% modelRef.nonlinearTolerance = 1.0e-9;
% modelRef.FacilityModel.toleranceWellBHP = 10;
% modelRef.FacilityModel.toleranceWellRate = 1.0e-9;
% modelRef.OutputStateFunctions{end+1} = 'ComponentTotalFlux';
% modelRef.OutputStateFunctions{end+1} = 'FaceComponentMobility';
% modelRef.OutputStateFunctions{end+1} = 'Mobility';
% 
% % 1.5 Run reference simulation
% [wellSolsRef, statesRef] = simulateScheduleAD(state0, modelRef, scheduleRef);

