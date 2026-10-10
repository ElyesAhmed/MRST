

%% Adjoint-Based Parameter Inversion for Reservoir Characterization
% This example demonstrates parameter inversion using adjoint optimization to
% reconstruct missing reservoir parameters by matching production data.
% The methodology follows Kohn-Vogelius optimization principles.

clear all;
close all;
clc;

%% Load Required MRST Modules
mrstModule add agglom upscaling coarsegrid ...
               ad-core ad-blackoil ad-props ...
               optimization deckformat
mrstModule add ad-core ad-blackoil ad-props optimization spe10
%% (1) Reference Exact Model Setup and Simulation
% This section creates and simulates the reference model with known parameters
% that we will try to reconstruct through inversion.

% 1.1 Create synthetic reference model with known properties
setupModel3DForReconstruction
%poro = smoothdata(modelRef.rock.poro);
% Create model-object of class TwoPhaseOilWaterModel
modelRef = GenericBlackOilModel(G.Matrix,G.Matrix.rock, fluid, 'gas', false);
% Set initial state and run simulation:
state0 = initResSol(G.Matrix, 100*barsa, [0, 1]);
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
scheduleRef = perturbedSimpleSchedule(time_steps, 'W', W, ...
    'pressureFac', 0.375, 'rateFac', 0.45, 'perturbStep', perturbStep);

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
[wellSolsRef, statesRef] = simulateScheduleAD(state0, modelRef, scheduleRef);
scheduleRefPert = scheduleRef;
% perturbFactor = 0.02;
% for i = 1:length(scheduleRefPert.control)
%     for j = 1:2
%         originalVal = scheduleRefPert.control(i).W(j).val;
%         scheduleRefPert.control(i).W(j).val = originalVal * (1 + perturbFactor);
%     end
% 
%        
%     for j = 3:4
%         originalVal = scheduleRefPert.control(i).W(j).val;
%         scheduleRefPert.control(i).W(j).val = originalVal * (1 - perturbFactor);
%     end
% end
% 1.6 Visualize reference results
subplot(1,2,2); cla;
plotCellData(modelRef.G, modelRef.rock.poro, 'EdgeColor', 'none');
title('Reference Pore Volume');
plotFaces(modelRef.G, boundaryFaces(modelRef.G), ...
    'EdgeColor', [0.4 0.4 0.4], 'EdgeAlpha', 0.5, 'FaceColor', 'none');
plotWell(modelRef.G, W, 'Color', 'k'); view(3); axis tight off;

%% (2) Duplication Procedure Setup
% This section establishes the dual control systems required for the
% Kohn-Vogelius optimization approach and verifies their equivalence.

% 2.1 Setup boundary conditions for inversion
bcfaces = boundaryFaceIndices(G,'Left');
%bcfaces = boundaryFaces(G);
% bcfaces = [bcfaces;boundaryFaceIndices(G,'RIGHT')];
% bcfaces = [bcfaces;boundaryFaceIndices(G,'South')];
% bcfaces = [bcfaces;boundaryFaceIndices(G,'North')];

data = setupReservoirBoundaryConditions(modelRef, statesRef, scheduleRef, 'seismic',bcfaces);


% 2.2 Create original BHP-Rate control system
scheduleBhpRate = scheduleRefPert;
modelBhpRate = modelRef;
[wellSolsBhpRate, statesBhpRate] = simulateScheduleAD(state0, modelBhpRate, scheduleBhpRate);

% 2.3 Create converted Rate-BHP control system
convMap = struct('bhp', 'rate', 'rate', 'bhp');
scheduleRateBhp = convertWellControls(scheduleRefPert, statesRef, modelRef, ...
                    'ConversionMap', convMap, 'bc', data.bc);
[wellSolsRateBhp, statesRateBhp] = simulateScheduleAD(state0, modelBhpRate, scheduleRateBhp);

% 2.4 Validate control system equivalence
tolerance = 1e-3;
isValid = true;

for i = 1:numel(wellSolsBhpRate)
    % Compare BHP where applicable
    for w = 1:numel(wellSolsBhpRate{i})
        if isfield(wellSolsBhpRate{i}(w), 'bhp') && isfield(wellSolsRateBhp{i}(w), 'bhp')
            bhpDiff = abs(wellSolsBhpRate{i}(w).bhp - wellSolsRateBhp{i}(w).bhp)/barsa;
            if bhpDiff > tolerance
                fprintf('BHP mismatch at step %d, well %d: %.2e\n', i, w, bhpDiff);
                isValid = false;
            end
        end
        
        % Compare rates where applicable
        if isfield(wellSolsBhpRate{i}(w), 'qWs') && isfield(wellSolsRateBhp{i}(w), 'qWs')
            rateDiff = abs(wellSolsBhpRate{i}(w).qWs - wellSolsRateBhp{i}(w).qWs);
            if rateDiff > tolerance
                fprintf('Water rate mismatch at step %d, well %d: %.2e\n', i, w, rateDiff);
                isValid = false;
            end
        end
    end
end

if isValid
    fprintf('\nValidation PASSED: Both control systems produce equivalent results within tolerance.\n');
% else
%     error('Validation FAILED: Control systems produce different results');
end

% Apply consistent 2% perturbation to both schedules
% perturbFactor = 0.02;
% for i = 1:length(scheduleBhpRate.control)
%     for j = 1:2
%         originalVal = scheduleBhpRate.control(i).W(j).val;
%         scheduleBhpRate.control(i).W(j).val = originalVal * (1 + perturbFactor);
%     end
% 
%        
%     for j = 3:4
%         originalVal = scheduleBhpRate.control(i).W(j).val;setupModel3DForReconstruction
%         scheduleBhpRate.control(i).W(j).val = originalVal * (1 - perturbFactor);
%     end
% end
% 
% for i = 1:length(scheduleRateBhp.control)
%     for j = 1:2
%         originalVal = scheduleRateBhp.control(i).W(j).val;
%         scheduleRateBhp.control(i).W(j).val = originalVal * (1 - perturbFactor);
%     end
% 
%        
%     for j = :2
%         originalVal = scheduleRateBhp.control(i).W(j).val;
%         scheduleRateBhp.control(i).W(j).val = originalVal * (1 - perturbFactor);
%     end
% end

% 2.5 Visualize control system comparison
summary_plots2 = plotWellSols({wellSolsRef, wellSolsBhpRate}, ...
                            {scheduleRef.step.val, scheduleBhpRate.step.val}, ...
                            'datasetnames', {'Reference Model', 'BHP-Rate Model'});

summary_plots3 = plotWellSols({wellSolsBhpRate, wellSolsRateBhp}, ...
                            {scheduleBhpRate.step.val, scheduleRateBhp.step.val}, ...
                            'datasetnames', {'BHP-Rate Model', 'Rate-BHP Model'});

%% (3) Initial Model Setup for Optimization
% This section prepares the initial guess model that will be used to start
% the optimization procedure.
% 3.0 We construct a prior field
prior.location = data.wellCells;
prior.name = 'porevolume';
prior.value = modelRef.operators.pv(data.wellCells);


% 3.1 Initialize model with estimated porosity
 InitRock = modelRef.rock;
 % Create correlated Gaussian random field
rng(42); % For reproducibility
correlation_length = [10, 10, 1]; % In grid units
perturbation = gaussianField(G.cartDims, ...       % Grid dimensions [nx,ny,nz]
                           [-0.85, 0.85], ...       % Target interval for values
                           [5, 5, 1], ...        % Kernel size [5x5x1]
                           0.65);                % Standard deviation

% Add perturbation to base porosity
%InitRock.poro = InitRock.poro .* (1 + perturbation(:)).*0+0.2;
% Perturb reference rock slightly for calibration target
rock_per = rock;
%rock_per.perm = rock.perm * 1.1;  % 10% increase in permeability
%rock_per.perm(:,1) = rock_per.perm(:,1) .* (1 + 0.05*randn(size(rock_per.perm(:,1))));  % additional random pert
rock_per.poro = rock.poro * 0.9;
rock_per.poro = rock_per.poro .* (1 + 0.25*randn(size(rock_per.poro)));
rock_per.poro = max(min(rock_per.poro, 0.4025), 0.038);  % physical bounds
%InitRock.poro = min(max(InitRock.poro, 0.03), 0.41); % Keep within physical bounds
% initialPoroValue = 0.223;
% InitRock.poro = ones(size(InitRock.poro)) * initialPoroValue;
%InitRock.poro = InitRock.poro.*0 +  6.01.*min(InitRock.perm);

% 3.2 Create initial model with guessed properties
InitModel = GenericBlackOilModel(modelRef.G, rock_per, modelRef.fluid, ...
    'gas', false, 'water', true, 'oil', true);
InitModel.OutputStateFunctions{end+1} = 'ComponentTotalFlux';
InitModel.OutputStateFunctions{end+1} = 'FaceComponentMobility';
InitModel.OutputStateFunctions{end+1} = 'Mobility';

% 3.3 Preserve exact pore volume in well cells from seismic data
seismicdata.cells = data.wellCells;
seismicdata.value = modelRef.operators.pv(seismicdata.cells);
InitModel.operators.pv(seismicdata.cells) = seismicdata.value;

% 3.4 Configure model solver settings for optimization
InitModel.useCNVConvergence = false;
InitModel.nonlinearTolerance = 1.0e-9;
InitModel.FacilityModel.toleranceWellBHP = 1000;
InitModel.FacilityModel.toleranceWellRate = 1.0e-8;
InitModel.FacilityModel.nonlinearTolerance = 1.0e-8;
InitModel.OutputStateFunctions{end+1} = 'ComponentTotalFlux';

% 3.5 Validate initial model with both control systems
modelBhpRate = InitModel;
modelRateBhp = InitModel;
% for i = 1:18
% scheduleBhpRate.control(i).W(1).WI = 1.e-12;
% scheduleBhpRate.control(i).W(2).WI = 1.e-12;
% scheduleBhpRate.control(i).W(3).WI = 1.e-12;   
% scheduleBhpRate.control(i).W(4).WI =  1.e-12;  
% end
% 
% for i = 1:37
% scheduleRateBhp.control(i).W(1).WI = 1.e-12;
% scheduleRateBhp.control(i).W(2).WI = 1.e-12;
% scheduleRateBhp.control(i).W(3).WI = 1.e-12;   
% scheduleRateBhp.control(i).W(4).WI =  1.e-12;  
% end
fprintf('\nRunning initial forward simulations with guessed parameters...\n');
[wellSolsBhpRate, statesBhpRate] = simulateScheduleAD(state0, modelBhpRate, scheduleBhpRate);
[wellSolsRateBhp, statesRateBhp] = simulateScheduleAD(state0, modelRateBhp, scheduleRateBhp);

% 3.6 Visualize initial model performance
figure(3);
summary_plots3 = plotWellSols({wellSolsBhpRate, wellSolsRateBhp}, ...
                            {scheduleBhpRate.step.val, scheduleRateBhp.step.val},...
                            'datasetnames', {'BHP-Rate Model (Initial)', ...
                            'Rate-BHP Model (Initial)'});
title('Initial Model Performance Comparison');

%% Optimization Setup Complete
fprintf('\nInitial model setup complete. Ready to begin optimization procedure.\n');
%% Optimization Setup
% Prepare parameter structures for both control systems
params_modelBhpRate = [];  
params_modelRateBhp = [];

% Define model setups
setup_modelBhpRate = struct('model', InitModel, 'schedule', scheduleBhpRate, 'state0', state0);
setup_modelRateBhp = struct('model', InitModel, 'schedule', scheduleRateBhp, 'state0', state0);

%Set parameter bounds based on reference model
pv_range = [min(modelRef.operators.pv), max(modelRef.operators.pv)];
%trans_range = ([min(modelRef.operators.T), max(modelRef.operators.T)]);
params_modelBhpRate = addParameter(params_modelBhpRate, setup_modelBhpRate, ...
                                 'name', 'porevolume', 'boxLims', pv_range);
params_modelRateBhp = addParameter(params_modelRateBhp, setup_modelRateBhp, ...
                                 'name', 'porevolume', 'boxLims', pv_range);
% % 
% % trans_range = [1.0e-13, 2.0e-10];
% trans_range = ([min(modelRef.operators.T).*0.95, max(modelRef.operators.T).*2.05]);
% params_modelBhpRate = addParameter(params_modelBhpRate, setup_modelBhpRate, ...
%                                  'name', 'transmissibility', 'scaling',   'log',      'boxLims', trans_range);
% params_modelRateBhp = addParameter(params_modelRateBhp, setup_modelRateBhp, ...
%                                  'name', 'transmissibility',  'scaling',  'log', 'boxLims', trans_range);

%% Optimization Configuration
% Define weighting for objective function
weighting = {'WaterRateWeight', (500*stb/day)^-1, ...
             'OilRateWeight',   (500*stb/day)^-1, ...
             'BHPWeight',       (100*barsa)^-1};


%% Two-Stage Optimization Process
% 1. BFGS Optimization (Optional)
use_bfgs = true;  % Set to false to skip BFGS stage
max_iterations = 30;
bfgs_iterations = 5;

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

% 2. Levenberg-Marquardt Optimization
obj_func_lm = @(m1,s1,sch1,m2,s2,sch2,sc,obj_prior, deriv,t,s1_in,s2_in) ...
    matchObservedTwinModelsOW(m1,s1,sch1,m2,s2,sch2,sc, obj_prior,...
    'computePartials',deriv,'tstep',t,'state1',s1_in,...
    'state2',s2_in,'from_states',false,weighting{:}, 'mismatchSum', false);
fprintf('\nStarting Levenberg-Marquardt optimization...\n');
obj_lm = @(p) evaluateMatchTwinModelsSummands(p, obj_func_lm, setup_modelBhpRate, ...
                                            setup_modelRateBhp, params_modelBhpRate, ...
                                            statesRef, prior);

[~, p_opt_final, ~] = unitBoxLM(p_opt, obj_lm, 'maxIt', max_iterations);

%% Final Simulation with Optimized Parameters
setup_final = updateSetupFromScaledParameters(setup_modelBhpRate, params_modelBhpRate, p_opt_final);
setup_final.schedule.step = scheduleRef.step; % Use full time horizon

fprintf('\nRunning final simulation with optimized parameters...\n');
[wellSols_opt, states_opt] = simulateScheduleAD(setup_final.state0, ...
                                              setup_final.model, ...
                                              setup_final.schedule);

%% Results Visualization
plotWellSols({wellSolsRef, wellSols_opt}, ...
            {scheduleRef.step.val, setup_final.schedule.step.val}, ...
            'datasetnames', {'Reference Model', 'Optimized Model'});
title('Optimization Results Comparison');

% Create custom "boxplot" using percentiles
ref_pv = modelRef.operators.pv;
init_pv = InitModel.operators.pv;
opt_pv = setup_final.model.operators.pv;
percentiles = [5 25 50 75 95];
ref_prct = prctile(ref_pv, percentiles);
init_prct = prctile(init_pv, percentiles);
opt_prct = prctile(opt_pv, percentiles);

figure;
hold on;
% Reference
plot([1 1], [ref_prct(1) ref_prct(5)], 'b-', 'LineWidth', 2);
rectangle('Position', [0.8, ref_prct(2), 0.4, ref_prct(4)-ref_prct(2)], ...
          'FaceColor', [0.7 0.7 1], 'EdgeColor', 'b');
plot([0.8 1.2], [ref_prct(3) ref_prct(3)], 'b-', 'LineWidth', 2);

% Initial
plot([2 2], [init_prct(1) init_prct(5)], 'r-', 'LineWidth', 2);
rectangle('Position', [1.8, init_prct(2), 0.4, init_prct(4)-init_prct(2)], ...
          'FaceColor', [1 0.7 0.7], 'EdgeColor', 'r');
plot([1.8 2.2], [init_prct(3) init_prct(3)], 'r-', 'LineWidth', 2);

% Optimized
plot([3 3], [opt_prct(1) opt_prct(5)], 'g-', 'LineWidth', 2);
rectangle('Position', [2.8, opt_prct(2), 0.4, opt_prct(4)-opt_prct(2)], ...
          'FaceColor', [0.7 1 0.7], 'EdgeColor', 'g');
plot([2.8 3.2], [opt_prct(3) opt_prct(3)], 'g-', 'LineWidth', 2);

set(gca, 'XTick', 1:3, 'XTickLabel', {'Reference', 'Initial', 'Optimized'});
ylabel('Pore Volume');
title('Distribution Comparison (5-25-50-75-95 percentiles)');
grid on;

% Panel 3: Spatial distribution
figure;
plot(modelRef.operators.pv, 'LineWidth', 2); hold on;
plot(InitModel.operators.pv, '--', 'LineWidth', 1.5);
plot(setup_final.model.operators.pv, 'LineWidth', 1.5);
legend('Reference', 'Initial Guess', 'Optimized', 'Location','best');
xlabel('Cell Index'); ylabel('Pore Volume');
title('Spatial Distribution');
grid on;

% Calculate basic statistics
ref_pv = modelRef.operators.pv;
init_pv = InitModel.operators.pv;
opt_pv = setup_final.model.operators.pv;

stats = struct();
stats.Reference = [mean(ref_pv), std(ref_pv), min(ref_pv), max(ref_pv)];
stats.Initial = [mean(init_pv), std(init_pv), min(init_pv), max(init_pv)];
stats.Optimized = [mean(opt_pv), std(opt_pv), min(opt_pv), max(opt_pv)];

% Display as formatted table
fprintf('%-12s %-8s %-8s %-8s %-8s\n', 'Model', 'Mean', 'Std', 'Min', 'Max');
fprintf('%-12s %-8.2f %-8.2f %-8.2f %-8.2f\n', 'Reference', stats.Reference);
fprintf('%-12s %-8.2f %-8.2f %-8.2f %-8.2f\n', 'Initial', stats.Initial);
fprintf('%-12s %-8.2f %-8.2f %-8.2f %-8.2f\n', 'Optimized', stats.Optimized);

% Plot comparison
figure;
bar([stats.Reference(1); stats.Initial(1); stats.Optimized(1)]);
set(gca, 'XTickLabel', {'Reference', 'Initial', 'Optimized'});
ylabel('Mean Pore Volume');
title('Mean PV Comparison');
grid on;