mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection reorder matlab_bgl

%%

gravity reset off;

n = 4;
l = 1000;
G = computeGeometry(cartGrid([n,n], [l,l]*meter));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [1, 1]                 );

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%

[jt, ot, mt] = deal(Inf);
jt = 0.2;
ot = 1e-3;
degree = 1;
disc   = DGDiscretization(modelDG.transportModel, 'degree', degree, ...
                         'basis', 'legendre', 'useUnstructCubature', true,  'jumpTolerance', jt, ...
     'outTolerance', ot, 'meanTolerance', mt);
modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);    

%%

time = 2*year;
rate = 1*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
state0 = assignDofFromState(modelDG.transportModel.disc, state0);

%%

[wsDG, statesDG, rep] = simulateScheduleAD(state0, modelDG, schedule);

%%
    
[modelDG.transportModel.extraStateOutput, modelDG.pressureModel.extraStateOutput] = deal(true);
modelDGreorder = modelDG;
modelDGreorder.pressureModel.extraStateOutput = true;

modelDGreorder.transportModel = ReorderingModelDG_ghost(modelDGreorder.transportModel, 'plotProgress', true);

modelDGreorder.transportModel.chunkSize = 1;
modelDGreorder.transportModel.parent.extraStateOutput = true;


[wsDGReorder, statesDGReorder, repDGReorder] = simulateScheduleAD(state0, modelDGreorder, schedule);

%%

close all

figure
plotToolbar(G, statesDG); colormap(jet);
figure


plotToolbar(G, statesDGReorder); colormap(jet);

plotWellSols({wsDG, wsDGReorder});

%%

mrstModule add upr
n = 3;
G = pebiGrid(l/n, [l,l]);
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);

W = [];
xw = [0,l/2; l,l/2];
isInj = [true, false];
for wNo = 1:size(xw,1)
    d = sqrt(sum((xw(wNo,:) - G.cells.centroids).^2, 2));
    wc = find(d == min(d));
    wc = wc(1);
    if isInj(wNo)
        W = addWell(W, G, rock, wc, 'type', 'rate', 'val', rate, 'comp_i', [1,0]);
    else
        W = addWell(W, G, rock, wc, 'type', 'bhp', 'val', 50*barsa, 'comp_i', [1,0]);
    end
end

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

%%
    
modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%

[jt, ot, mt] = deal(Inf);
% jt = 0.2;
% ot = 1e-3;
degree = 0;
disc   = DGDiscretization(modelDG.transportModel, 2, 'degree', degree, ...
                         'basis', 'legendre', 'useUnstructCubature', true,  'jumpTolerance', jt, ...
     'outTolerance', ot, 'meanTolerance', mt);
modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);


sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
state0 = assignDofFromState(modelDG.transportModel.disc, state0);

%%

[modelDG.transportModel.extraStateOutput, modelDG.pressureModel.extraStateOutput] = deal(true);
modelDGreorder = modelDG;
modelDGreorder.pressureModel.extraStateOutput = true;

modelDGreorder.transportModel = ReorderingModelDG_ghost(modelDGreorder.transportModel, 'plotProgress', false);

modelDGreorder.transportModel.chunkSize = 1;
modelDGreorder.transportModel.parent.extraStateOutput = true;


[wsDGReorder, statesDGReorder, repDGReorder] = simulateScheduleAD(state0, modelDGreorder, schedule);

%%

[wsDG, statesDG, rep] = simulateScheduleAD(state0, modelDG, schedule);