% Extract grid and rock from Eclipse deck
deck = readEclipseDeck('/home/elyes/Downloads/spe11-decks-main/csp11b/isothermal/130_62/CSP11B_DISGAS.DATA');
deck = convertDeckUnits(deck);
[state0,model,schedule,~] = initEclipseProblemAD(deck, 'getSchedule',true, 'getInitialState',true);

% Use the grid and rock from the Eclipse deck
G = model.G;
% [state0,model,schedule,~] = initEclipseProblemAD(deck, 'getSchedule',true, 'getInitialState',true);

rock = model.rock;

% fluid
pRef = 200*barsa;

fluid = initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
c = 0.*1e-5/barsa;
%pRef = 200*barsa;
fluid.bO = @(p) exp((p - pRef)*c);

W = [];
% Injectors (lower-left and upper-right)
ci(1) = schedule.control(3).W(1).cells;
ci(2) = schedule.control(3).W(2).cells;
for k  = 1:1
    W = addWell(W, G, rock, ci(k), 'Type' , 'rate', ...
                                   'Val'  , 2.*schedule.control(1).W(1).val, ...
                                   'Name' , sprintf('I%d', k), ...
                                   'comp_i', [1 0], ...
                                   'Sign' , 1);
end
% Producers (upper-left and -right)
cp(1) = G.cartDims(1);
cp(2) = 1 + (G.cartDims(2)-1)*G.cartDims(1);
for k  = 2:2
    W = addWell(W, G, rock, ci(k), 'Type', 'bhp', ...
                                   'Val' , 150*barsa, ...
                                   'Name', sprintf('P%d', k), ...
                                   'comp_i', [0 1], ...
                                   'Sign', -1);
end

% % Set up 4 control-steps each 150 days
ts = { [1 1 3 5 5 10 10 10 15 15 15 15 15 15 15]'*year, ...
                   repmat(150/10, 10, 1)*year, ...
                   repmat(150/6, 6, 1)*year, ...
                   repmat(150/6, 6, 1)*year};
       
numCnt = numel(ts);
[schedule.control(1:numCnt).W] = deal(W);
schedule.step.control = rldecode((1:4)', cellfun(@numel, ts));
schedule.step.val     = vertcat(ts{:});

gravity off

%poro = smoothdata(modelRef.rock.poro);
% Create model-object of class TwoPhaseOilWaterModel
modelRef = GenericBlackOilModel(G, rock, model.fluid, 'gas', false);
% Set initial state and run simulation:
%state0 = initResSol(G, state0.pressure, state0.s);
% 1.2 Visualize reference model
figure(1);
subplot(1,2,1);
plotCellData(modelRef.G, modelRef.rock.poro, 'EdgeAlpha', 0.5); 
view(3); axis tight off;
title('Reference Model Porosity');
plotWell(modelRef.G, W, 'Color', 'k');
drawnow;

% 1.3 Create and perturb schedule for realistic inversion scenario
time_steps = schedule.step.val;
W = schedule.control(1).W;

% Add controlled perturbations to well operations
rng(0); % For reproducibility
perturbStep = [ones(1,4), round(.5+(2:(numel(schedule.step.val)-3))/2)]';
scheduleRef = schedule;perturbedSimpleSchedule(time_steps, 'W', W, ...
    'pressureFac', 0.75, 'rateFac', 0.9, 'perturbStep', perturbStep);

% 1.4 Configure model for precise simulations
modelRef.useCNVConvergence = false; 
modelRef.toleranceMB = 1.0e-9;
modelRef.nonlinearTolerance = 1.0e-9;
modelRef.FacilityModel.toleranceWellBHP = 10;
modelRef.FacilityModel.toleranceWellRate = 1.0e-9;
modelRef.OutputStateFunctions{end+1} = 'ComponentTotalFlux';
modelRef.OutputStateFunctions{end+1} = 'FaceComponentMobility';
modelRef.OutputStateFunctions{end+1} = 'Mobility';

% 1.5 Run reference simulation
[wellSolsRef, statesRef] = simulateScheduleAD(state0, model, scheduleRef);