mrstModule add spe10 vem vemmech dg ad-core ad-props ad-blackoil ...
    blackoil-sequential incomp msrsb matlab_bgl coarsegrid mrst-gui ...
    reorder weno vista compositional

mrstVerbose on

gravity reset off

%% Common stuff

name = 'qfs_co2_small';
baseName = 'qfs-CO2';
dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);
[state0, model, schedule] = getReorderingCase(name);

pack = @(state0, model, name, desc) ...
    packSimulationProblem(state0, model, schedule, baseName, ...
                          'Name'           , name          , ...
                          'Directory'      , dataDir       , ...
                          'Description'    , desc          );


%% Fully-implicit

model.AutoDiffBackend = DiagonalAutoDiffBackend();

fim = pack(state0, model, 'fim', 'Fully-implicit');

%% Sequential-implicit

modelseq = getSequential_local(model);
modelseq.transportNonLinearSolver.LinearSolver = BackslashSolverAD();
seq = pack(state0, modelseq, 'seq', 'Sequential');

%% Sequential-reordering

modelreorder = getSequential_local(model);
modelreorder.transportNonLinearSolver.LinearSolver = BackslashSolverAD();
modelreorder.pressureModel.extraStateOutput = true;
modelreorder.transportModel = ReorderingModel(modelreorder.transportModel);
modelreorder.transportModel.parent.extraStateOutput = true;

modelreorder.transportModel.chunkSize = 10;
modelreorder.transportModel.buffer = 0;

reorder = pack(state0, modelreorder, 'reorder', 'Reordering');

%% Sequential-reordering strict tolerance

modelreorder = getSequential_local(model);
modelreorder.pressureModel.extraStateOutput = true;
modelreorder.transportModel = ReorderingModel(modelreorder.transportModel);
modelreorder.transportModel.parent.extraStateOutput = true;

modelreorder.transportModel.chunkSize = 10;
modelreorder.transportModel.buffer = 0;
modelreorder.transportNonLinearSolver.LinearSolver = BackslashSolverAD();
modelreorder.transportModel.parent.nonlinearTolerance = 1e-5;

reorderStrict = pack(state0, modelreorder, 'reorder-strict', 'Reordering Strict');

%% Adaptive

nls = NonLinearSolver();

modeladapt = getSequential_local(model);
modeladapt.pressureModel.extraStateOutput = true;

p = partitionCartGrid(model.G.cartDims, [15,15,1]);
G = model.G;
G = computeCellDimensions2(G);
GC = generateCoarseGrid(G, p);
m = getRefinementMappings(GC, GC, model.G, [1, GC.cells.num]);

GC = generateCoarseGrid(model.G, m.newPartition);
GC.cells.refined = m.refined;

coarsemodel = upscaleModelTPFA(modeladapt.transportModel, m.newPartition);
state0C = upscaleState(coarsemodel, modeladapt.transportModel, state0);
modeladapt = AdaptiveSequentialPressureTransportModel(modeladapt.pressureModel, modeladapt.transportModel, GC);
state0 = ensureDensitiesPresentCompositional(model, state0);
state0C.rho = repmat(state0.rho(1,:), GC.cells.num,1);
% state0C.bfactor = state0C.rho./[fluid.rhoOS, fluid.rhoGS];
% state0C.x = repmat(state0.x(1,:), GC.cells.num, 1);
% state0C.y = repmat(state0.y(1,:), GC.cells.num, 1);
% state0C.L = repmat(state0.L(1,:), GC.cells.num, 1);
% state0C.K = repmat(state0.K(1,:), GC.cells.num, 1);
% state0C.Z_V = repmat(state0.Z_V(1,:), GC.cells.num, 1);
% state0C.Z_L = repmat(state0.Z_L(1,:), GC.cells.num, 1);

state0C.components = repmat(state0C.components(1,:), GC.cells.num, 1);
state0.transportState = state0C;
state0.transportState.G = GC;
state0.transportState.pv = coarsemodel.operators.pv;
state0.transportModel = modeladapt.transportModel;
modeladapt.plotProgress = true;
modeladapt.refineTol = 1e-2;
modeladapt.coarsenTol = 0.5*1e-2;

adapt = pack(state0, modeladapt, 'adapt', 'Adaptive');

%%

coarsemodel = upscaleModelTPFA(model, m.newPartition);

fluid = model.fluid;
state0C.bfactor = state0C.rho./[fluid.rhoOS, fluid.rhoGS];
state0C.x = repmat(state0.x(1,:), GC.cells.num, 1);
state0C.y = repmat(state0.y(1,:), GC.cells.num, 1);
state0C.L = repmat(state0.L(1,:), GC.cells.num, 1);
state0C.K = repmat(state0.K(1,:), GC.cells.num, 1);
state0C.Z_V = repmat(state0.Z_V(1,:), GC.cells.num, 1);
state0C.Z_L = repmat(state0.Z_L(1,:), GC.cells.num, 1);

coarsemodel = getSequential_local(coarsemodel);

W = schedule.control(1).W;
WC = W;
d = pdist2(coarsemodel.transportModel.G.cells.centroids, model.G.cells.centroids([W.cells],:));
[~, c] = min(d);
for wNo = 1:numel(W)
    WC(wNo).cells = c(wNo);
end

coarse = pack(state0C, coarsemodel, 'coarse', 'Coarse');
coarse.SimulatorSetup.schedule.control(1).W = WC;

%%

problems = {fim, seq, reorder, adapt, coarse, reorderStrict};

%% Simulate problems

runIx = 2:4;
for pNo = runIx
    [ok, status] = simulatePackedProblem(problems{pNo});
end

%%

setup = problems{3}.SimulatorSetup;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);








