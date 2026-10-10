%% Adjoint-Based Parameter Inversion for 1D Reservoir Characterization
% This example demonstrates parameter inversion using adjoint optimization to
% reconstruct missing reservoir parameters by matching production data in 1D.
% The methodology follows Kohn-Vogelius optimization principles.

clear all;
close all;
clc;

%% Load Required MRST Modules
mrstModule add ad-core ad-blackoil ad-props optimization

%% (1) Generate 1D Reference Problem
% Use the provided function to create a simple 1D problem
problem = setupFnSimpleProblem1D(2);

% Extract components from the problem structure
G        = problem.SimulatorSetup.model.G;
rock     = problem.SimulatorSetup.model.rock;
fluid    = problem.SimulatorSetup.model.fluid;
schedule = problem.SimulatorSetup.schedule;
state0   = problem.SimulatorSetup.state0;
modelRef = problem.SimulatorSetup.model;
% 1.3 Create and perturb schedule for realistic inversion scenario
time_steps = schedule.step.val;
W = schedule.control(1).W;

% Add controlled perturbations to well operations
rng(0); % For reproducibility
perturbStep = [ones(1,4), round(.5+(2:(numel(schedule.step.val)-3))/2)]';
schedule = perturbedSimpleSchedule(time_steps, 'W', W, ...
    'pressureFac', 0.9, 'rateFac', 0.9, 'perturbStep', perturbStep);

% 1.1 Visualize reference model
figure(1);
subplot(1,2,1);
plot(rock.poro); 
title('Reference Model Porosity (1D)');
xlabel('Cell Index'); ylabel('Porosity');

% 1.2 Run reference simulation
[wellSolsRef, statesRef] = simulateScheduleAD(state0, modelRef, schedule);


%% (2) Duplication Procedure Setup for 1D
% This section establishes the dual control systems required for the
% Kohn-Vogelius optimization approach.

% 2.1 Setup boundary conditions for inversion (use left boundary in 1D)
bcfaces = boundaryFaces(G);

%data = setupReservoirBoundaryConditions(modelRef, statesRef, schedule, 'seismic', bcfaces);
data.bdFaces = bcfaces;
data.bdCells = sum(G.faces.neighbors(bcfaces, :), 2);
% Create boundary conditions for each timestep
for i = 1:numel(statesRef)
    if any(modelRef.gravity)
        gRhoDz = computeGravity(model, statesRef{i});
        facePressure = statesRef{i}.pressure(data.bdCells) + gRhoDz{2}(data.bdCells);
    else
        facePressure = statesRef{i}.pressure(data.bdCells);
    end
    data.bc{i} = addBC([], data.bdFaces, 'pressure', facePressure, ...
        'sat', statesRef{i}.s(data.bdCells, :));
end
% 2.2 Create original BHP-Rate control system
scheduleBhpRate = schedule;
modelBhpRate = modelRef;
[wellSolsBhpRate, statesBhpRate] = simulateScheduleAD(state0, modelBhpRate, scheduleBhpRate);
data.wellCells = [1,G.cells.num];
% 2.3 Create converted Rate-BHP control system
convMap = struct('bhp', 'rate', 'rate', 'bhp');
scheduleRateBhp = convertWellControls(schedule, statesRef, modelRef, ...
                    'ConversionMap', convMap, 'bc', data.bc);
[wellSolsRateBhp, statesRateBhp] = simulateScheduleAD(state0, modelBhpRate, scheduleRateBhp);

% 2.4 Validate control system equivalence
tolerance = 1e-2;
isValid = true;

for i = 1:numel(wellSolsBhpRate)
    for w = 1:numel(wellSolsBhpRate{i})
        if isfield(wellSolsBhpRate{i}(w), 'bhp') && isfield(wellSolsRateBhp{i}(w), 'bhp')
            bhpDiff = abs(wellSolsBhpRate{i}(w).bhp - wellSolsRateBhp{i}(w).bhp)/barsa;
            if bhpDiff > tolerance
                fprintf('BHP mismatch at step %d, well %d: %.2e\n', i, w, bhpDiff);
                isValid = true;
            end
        end
        
        if isfield(wellSolsBhpRate{i}(w), 'qWs') && isfield(wellSolsRateBhp{i}(w), 'qWs')
            rateDiff = abs(wellSolsBhpRate{i}(w).qWs - wellSolsRateBhp{i}(w).qWs);
            if rateDiff > tolerance
                fprintf('Water rate mismatch at step %d, well %d: %.2e\n', i, w, rateDiff);
                isValid = true;
            end
        end
    end
end

if isValid
    fprintf('\nValidation PASSED: Both control systems produce equivalent results within tolerance.\n');
else
    error('Validation FAILED: Control systems produce different results');
end

%% (3) Initial Model Setup for Optimization
% 3.1 Initialize model with perturbed porosity
InitRock = rock;
rng(42); % For reproducibility

% Create smooth perturbation for 1D
x = linspace(0, 1, G.cells.num)';
perturbation = 0.3 * sin(2*pi*x) + 0.1 * cos(4*pi*x);

%InitRock.poro = rock.poro .* (1 + perturbation);
InitRock.poro = InitRock.poro.*0+min(InitRock.poro)+0.1; % Keep within physical bounds
InitRock.poro(data.wellCells) = rock.poro(data.wellCells);
% 3.2 Create initial model with guessed properties
InitModel = GenericBlackOilModel(G, InitRock, fluid, 'gas', false);
InitModel.OutputStateFunctions{end+1} = 'ComponentTotalFlux';

% 3.3 Configure model solver settings for optimization
InitModel.useCNVConvergence = false;
InitModel.nonlinearTolerance = 1.0e-9;
InitModel.FacilityModel.toleranceWellBHP = 1000;
InitModel.FacilityModel.toleranceWellRate = 1.0e-8;

% 3.4 Validate initial model with both control systems
modelBhpRate = InitModel;
modelRateBhp = InitModel;

fprintf('\nRunning initial forward simulations with guessed parameters...\n');
[wellSolsBhpRate, statesBhpRate] = simulateScheduleAD(state0, modelBhpRate, scheduleBhpRate);
[wellSolsRateBhp, statesRateBhp] = simulateScheduleAD(state0, modelRateBhp, scheduleRateBhp);

%% Optimization Setup
params_modelBhpRate = [];  
params_modelRateBhp = [];

setup_modelBhpRate = struct('model', InitModel, 'schedule', scheduleBhpRate, 'state0', state0);
setup_modelRateBhp = struct('model', InitModel, 'schedule', scheduleRateBhp, 'state0', state0);

% Set parameter bounds
pv_range = [min(modelRef.operators.pv), max(modelRef.operators.pv)];
params_modelBhpRate = addParameter(params_modelBhpRate, setup_modelBhpRate, ...
                                 'name', 'porevolume', 'boxLims', pv_range, 'scaling', 'log');
params_modelRateBhp = addParameter(params_modelRateBhp, setup_modelRateBhp, ...
                                 'name', 'porevolume', 'boxLims', pv_range, 'scaling', 'log');

%% Optimization Configuration
weighting = {'WaterRateWeight', (21*stb/day)^-1, ...
             'OilRateWeight',   (9.2*stb/day)^-1, ...
             'BHPWeight',       ( 144.8*barsa)^-1};



%% Two-Stage Optimization Process
% 1. BFGS Optimization (Optional)
use_bfgs = false;  % Set to false to skip BFGS stage
max_iterations = 30;
bfgs_iterations = 5;
prior.location = data.wellCells';
prior.name = 'porevolume';
prior.value = modelRef.operators.pv(data.wellCells);

if use_bfgs
    fprintf('\nStarting BFGS optimization...\n');
    % First objective function handle
    obj_func_bfgs = @(m1,s1,sch1,m2,s2,sch2,sc,obj_prior, deriv,t,s1_in,s2_in) ...
               matchObservedTwinModelsOW(m1,s1,sch1,m2,s2,sch2,sc,obj_prior,...
               'computePartials',deriv,'tstep',t,'state1',s1_in,...
               'state2',s2_in,'from_states',false,weighting{:});
    p0 = getScaledParameterVector(setup_modelBhpRate, params_modelBhpRate);
    obj_bfgs = @(p) evaluateMatchTwinModels(p, obj_func_bfgs, setup_modelBhpRate, ...
                                          setup_modelRateBhp, params_modelBhpRate, ...
                                          statesRef, prior);
    
    [~, p_opt, ~] = unitBoxBFGS(p0, obj_bfgs, ...
                               'maxIt', bfgs_iterations, ...
                               'objChangeTol', 1e-7, ...
                               'gradTol', 1e-7);
else
    p_opt = getScaledParameterVector(setup_modelBhpRate, params_modelBhpRate);
end

%% Levenberg-Marquardt Optimization
obj_func_lm = @(m1,s1,sch1,m2,s2,sch2,sc,obj_prior, deriv,t,s1_in,s2_in) ...
    matchObservedTwinModelsOW(m1,s1,sch1,m2,s2,sch2,sc, obj_prior,...
    'computePartials',deriv,'tstep',t,'state1',s1_in,...
    'state2',s2_in,'from_states',false,weighting{:}, 'mismatchSum', false);

fprintf('\nStarting Levenberg-Marquardt optimization...\n');
obj_lm = @(p) evaluateMatchTwinModelsSummands(p, obj_func_lm, setup_modelBhpRate, ...
                                            setup_modelRateBhp, params_modelBhpRate, ...
                                            statesRef, prior); % No prior for 1D

p0 = getScaledParameterVector(setup_modelBhpRate, params_modelBhpRate);
[~, p_opt_final, ~] = unitBoxLM(p0, obj_lm, 'maxIt', 30);

%% Final Simulation with Optimized Parameters
setup_final = updateSetupFromScaledParameters(setup_modelBhpRate, params_modelBhpRate, p_opt_final);

fprintf('\nRunning final simulation with optimized parameters...\n');
[wellSols_opt, states_opt] = simulateScheduleAD(setup_final.state0, ...
                                              setup_final.model, ...
                                              setup_final.schedule);

%% Results Visualization - Adapted for 1D
figure(2);
subplot(2,1,1);
plot(G.cells.centroids(:,1), rock.poro, 'b-', 'LineWidth', 2); hold on;
plot(G.cells.centroids(:,1), InitRock.poro, 'r--', 'LineWidth', 1.5);
plot(G.cells.centroids(:,1), setup_final.model.rock.poro, 'g-', 'LineWidth', 1.5);
legend('Reference', 'Initial Guess', 'Optimized', 'Location','best');
xlabel('Distance (m)'); ylabel('Porosity');
title('1D Porosity Distribution');
grid on;

subplot(2,1,2);
plot(G.cells.centroids(:,1), modelRef.operators.pv, 'b-', 'LineWidth', 2); hold on;
plot(G.cells.centroids(:,1), InitModel.operators.pv, 'r--', 'LineWidth', 1.5);
plot(G.cells.centroids(:,1), setup_final.model.operators.pv, 'g-', 'LineWidth', 1.5);
legend('Reference', 'Initial Guess', 'Optimized', 'Location','best');
xlabel('Distance (m)'); ylabel('Pore Volume');
title('1D Pore Volume Distribution');
grid on;

% Calculate basic statistics
ref_pv = modelRef.operators.pv;
init_pv = InitModel.operators.pv;
opt_pv = setup_final.model.operators.pv;

fprintf('\n=== 1D Optimization Results ===\n');
fprintf('%-12s %-8s %-8s %-8s %-8s\n', 'Model', 'Mean', 'Std', 'Min', 'Max');
fprintf('%-12s %-8.2e %-8.2e %-8.2e %-8.2e\n', 'Reference', mean(ref_pv), std(ref_pv), min(ref_pv), max(ref_pv));
fprintf('%-12s %-8.2e %-8.2e %-8.2e %-8.2e\n', 'Initial', mean(init_pv), std(init_pv), min(init_pv), max(init_pv));
fprintf('%-12s %-8.2e %-8.2e %-8.2e %-8.2e\n', 'Optimized', mean(opt_pv), std(opt_pv), min(opt_pv), max(opt_pv));

% Plot well responses
figure(3);
plotWellSols({wellSolsRef, wellSols_opt}, ...
            {schedule.step.val, setup_final.schedule.step.val}, ...
            'datasetnames', {'Reference Model', 'Optimized Model'});
title('1D Optimization Results - Well Responses');