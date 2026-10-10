%% Adjoint-Based Parameter Inversion with Noise on Injection Wells
% Goal: Reconstruct missing production well BHP (wells 3 & 4)
% using noisy observations from injection wells (wells 1 & 2).
clear all; close all; clc;

%% Load Required MRST Modules
mrstModule add agglom upscaling coarsegrid ...
               ad-core ad-blackoil ad-props ...
               optimization deckformat spe10

%% (1) Reference Model Setup
setupModel2DForReconstruction
modelRef = GenericBlackOilModel(G, rock, fluid, 'gas', false);
state0 = initResSol(G, 100*barsa, [0, 1]);

% Schedule with perturbations (original true controls)
time_steps = schedule.step.val;
W = schedule.control(1).W;
rng(0);
perturbStep = [ones(1,4), round(.5+(2:(numel(schedule.step.val)-3))/2)]';
scheduleRef = perturbedSimpleSchedule(time_steps, 'W', W, ...
    'pressureFac', 0.375, 'rateFac', 0.45, 'perturbStep', perturbStep);

% Model configuration
modelRef.useCNVConvergence = false;
modelRef.toleranceMB = 1.0e-9;
modelRef.nonlinearTolerance = 1.0e-9;
modelRef.FacilityModel.toleranceWellBHP = 10;
modelRef.FacilityModel.toleranceWellRate = 1.0e-9;
modelRef.OutputStateFunctions{end+1} = 'ComponentTotalFlux';
modelRef.OutputStateFunctions{end+1} = 'FaceComponentMobility';
modelRef.OutputStateFunctions{end+1} = 'Mobility';

% Reference simulation (true model, true controls)
[wellSolsTrue, statesTrue] = simulateScheduleAD(state0, modelRef, scheduleRef);

% Visualisation
figure(1);
subplot(1,2,1);
plotCellData(modelRef.G, modelRef.rock.poro, 'EdgeAlpha', 0.5); view(3); axis tight off;
title('Reference Porosity'); plotWell(modelRef.G, W, 'Color', 'k');
subplot(1,2,2);
plotCellData(modelRef.G, modelRef.rock.poro, 'EdgeColor', 'none');
title('Reference Pore Volume'); plotFaces(modelRef.G, boundaryFaces(modelRef.G), ...
    'EdgeColor', [0.4 0.4 0.4], 'EdgeAlpha', 0.5, 'FaceColor', 'none');
plotWell(modelRef.G, W, 'Color', 'k'); view(3); axis tight off;

%% (2) Dual Control Systems Setup (original, noise‑free schedules)
bcfaces = boundaryFaceIndices(G,'Left');
data = setupReservoirBoundaryConditions(modelRef, statesTrue, scheduleRef, 'wells');

scheduleBhpRate_base = scheduleRef;
scheduleRateBhp_base = scheduleRef;
% Convert to Rate-BHP system (just for compatibility)
convMap = struct('bhp', 'bhp', 'rate', 'rate');
scheduleBhpRate_base = convertWellControls(scheduleRef, statesTrue, modelRef, ...
                    'ConversionMap', convMap, 'bc', []);
convMap = struct('bhp', 'bhp', 'rate', 'bhp');
scheduleRateBhp_base = convertWellControls(scheduleRef, statesTrue, modelRef, ...
                    'ConversionMap', convMap, 'bc', data.bc);

% Store true production well BHP (to be recovered)
trueBhpP1 = zeros(numel(scheduleRateBhp_base.control),1);
trueBhpP2 = zeros(numel(scheduleRateBhp_base.control),1);
for i = 1:numel(scheduleRateBhp_base.control)
    trueBhpP1(i) = scheduleRateBhp_base.control(i).W(3).val / barsa;
    trueBhpP2(i) = scheduleRateBhp_base.control(i).W(4).val / barsa;
end

% True water/oil rates for production wells (for plotting)
trueWaterRateP1 = zeros(numel(wellSolsTrue),1);
trueWaterRateP2 = zeros(numel(wellSolsTrue),1);
trueOilRateP1   = zeros(numel(wellSolsTrue),1);
trueOilRateP2   = zeros(numel(wellSolsTrue),1);
for t = 1:numel(wellSolsTrue)
    trueWaterRateP1(t) = abs(wellSolsTrue{t}(3).qWs);
    trueWaterRateP2(t) = abs(wellSolsTrue{t}(4).qWs);
    trueOilRateP1(t)   = abs(wellSolsTrue{t}(3).qOs);
    trueOilRateP2(t)   = abs(wellSolsTrue{t}(4).qOs);
end

%% (3) Initial Model Setup (rock properties guess)
InitRock = modelRef.rock;
rng(42);
InitModel = GenericBlackOilModel(modelRef.G, InitRock, modelRef.fluid, ...
    'gas', false, 'water', true, 'oil', true);
InitModel.OutputStateFunctions{end+1} = 'ComponentTotalFlux';
InitModel.OutputStateFunctions{end+1} = 'FaceComponentMobility';
InitModel.OutputStateFunctions{end+1} = 'Mobility';
seismicdata.cells = data.wellCells;
seismicdata.value = modelRef.operators.pv(seismicdata.cells);
InitModel.operators.pv(seismicdata.cells) = seismicdata.value;

InitModel.useCNVConvergence = false;
InitModel.nonlinearTolerance = 1.0e-9;
InitModel.FacilityModel.toleranceWellBHP = 1000;
InitModel.FacilityModel.toleranceWellRate = 1.0e-8;
InitModel.FacilityModel.nonlinearTolerance = 1.0e-8;

% Initial guess for production well BHP (set to 100 bar)
% We will use a copy of the base schedule and modify only production wells
scheduleBhpRate_init = scheduleBhpRate_base;
scheduleRateBhp_init = scheduleRateBhp_base;
for i = 1:length(scheduleBhpRate_init.control)
    scheduleBhpRate_init.control(i).W(3).val = 100*barsa();
    scheduleBhpRate_init.control(i).W(4).val = 100*barsa();
end
for i = 1:length(scheduleRateBhp_init.control)
    scheduleRateBhp_init.control(i).W(3).val = 100*barsa();
    scheduleRateBhp_init.control(i).W(4).val = 100*barsa();
end

% Initial forward runs with guessed parameters (using the base injection controls)
[wellSolsBhpRate_init, ~] = simulateScheduleAD(state0, InitModel, scheduleBhpRate_init);
[wellSolsRateBhp_init, ~] = simulateScheduleAD(state0, InitModel, scheduleRateBhp_init);
figure(2);
plotWellSols({wellSolsBhpRate_init, wellSolsRateBhp_init}, ...
             {scheduleBhpRate_init.step.val, scheduleRateBhp_init.step.val},...
             'datasetnames', {'BHP-Rate Model (Initial)', 'Rate-BHP Model (Initial)'});
title('Initial Model Performance Comparison');

%% Define prior field (for regularisation)
prior.location = data.wellCells;
prior.name = 'porevolume';
prior.value = modelRef.operators.pv(data.wellCells);

%% Weighting for objective
weighting = {'WaterRateWeight', (300*stb/day)^-1, ...
             'OilRateWeight',   (300*stb/day)^-1, ...
             'BHPWeight',       (150*barsa)^-1};

%% (4) Noise Levels to Test (applied to INJECTION wells, wells 1 & 2)
noiseLevels = [0, 0.01, 0.02, 0.05, 0.1];
numNoise = length(noiseLevels);

% Preallocate storage
objHistories = cell(numNoise,1);
finalModels = cell(numNoise,1);
finalParams = cell(numNoise,1);
recoveredBhpP1 = cell(numNoise,1);
recoveredBhpP2 = cell(numNoise,1);
recoveredWaterRateP1 = cell(numNoise,1);
recoveredWaterRateP2 = cell(numNoise,1);
recoveredOilRateP1 = cell(numNoise,1);
recoveredOilRateP2 = cell(numNoise,1);

%% (5) Loop Over Noise Levels
for iNoise = 1:numNoise
    relNoise = noiseLevels(iNoise);
    fprintf('\n========================================================\n');
    fprintf('TESTING NOISE LEVEL: %.1f%% (applied to INJECTION wells BHP)\n', relNoise*100);
    fprintf('========================================================\n');
    
    % ---- Create perturbed schedules for injection wells (1 & 2) ----
    scheduleBhpRate_pert = scheduleBhpRate_init;
    scheduleRateBhp_pert = scheduleRateBhp_init;
    rng(0); % reproducible noise
    for iCtrl = 1:numel(scheduleBhpRate_pert.control)
        for wIdx = 1:2   % wells 1 and 2 (injection)
            trueVal_BhpRate = scheduleBhpRate_init.control(iCtrl).W(wIdx).val;
            noise = trueVal_BhpRate * relNoise;
            scheduleBhpRate_pert.control(iCtrl).W(wIdx).val = trueVal_BhpRate + noise;
            %trueVal_RateBhp = scheduleRateBhp_base.control(iCtrl).W(wIdx).val;
            %noise = trueVal_RateBhp * relNoise * randn();
            %scheduleRateBhp_pert.control(iCtrl).W(wIdx).val = trueVal_RateBhp + noise;
        end
        % Production wells (3,4) keep true values for generating observations
    end

        
    for iCtrl = 1:numel(scheduleRateBhp_pert.control)    
        for wIdx = 1:2   % wells 1 and 2 (injection)
            %trueVal_BhpRate = scheduleBhpRate_base.control(iCtrl).W(wIdx).val;
            %noise = trueVal_BhpRate * relNoise * randn();
            %scheduleBhpRate_pert.control(iCtrl).W(wIdx).val = trueVal_BhpRate + noise;
            trueVal_RateBhp = scheduleRateBhp_init.control(iCtrl).W(wIdx).val;
            noise = trueVal_RateBhp * relNoise;
            scheduleRateBhp_pert.control(iCtrl).W(wIdx).val = trueVal_RateBhp + noise;
        end
        % Production wells (3,4) keep true values for generating observations
    end
    
    % ---- Generate observed data by simulating perturbed schedules ----
    [~, statesObserved] = simulateScheduleAD(state0, modelRef, scheduleBhpRate_pert);
    
    % ---- Create optimization setups with perturbed injection controls ----
    % We keep the same initial guess for production wells (100 bar)
    setup_BhpRate = struct('model', InitModel, 'schedule', scheduleBhpRate_pert, 'state0', state0);
    setup_RateBhp = struct('model', InitModel, 'schedule', scheduleRateBhp_pert, 'state0', state0);
    
    % However, for the forward simulation during optimization, we must use the
    % perturbed injection well controls. So we replace the injection well values
    % in the schedules of the setups with the perturbed ones.
    % The production well values will be overwritten by the optimizer.
%     for iCtrl = 1:numel(setup_BhpRate.schedule.control)
%         for wIdx = 1:2
%             setup_BhpRate.schedule.control(iCtrl).W(wIdx).val = scheduleBhpRate_pert.control(iCtrl).W(wIdx).val;
%             setup_RateBhp.schedule.control(iCtrl).W(wIdx).val = scheduleRateBhp_pert.control(iCtrl).W(wIdx).val;
%         end
%     end
    
    % ---- Define parameters (production well BHP) for this setup ----
    params_BhpRate = [];
    params_RateBhp = [];
    for i=1:37
        params_BhpRate = addParameter(params_BhpRate, setup_BhpRate, ...
            'name', 'bhp', 'type', 'value', 'belongsTo', 'schedule', ...
            'location', {'control', 'W', 'val'}, 'subset', [3,4], ...
            'controlType', 'bhp', 'boxLims', [80*barsa, 220*barsa], ...
            'scaling', 'linear', 'controlSteps', i);
        params_RateBhp = addParameter(params_RateBhp, setup_RateBhp, ...
            'name', 'rate', 'type', 'value', 'belongsTo', 'schedule', ...
            'location', {'control', 'W', 'val'}, 'subset', [3,4], ...
            'boxLims', [80*barsa, 220*barsa], 'controlType', 'bhp', ...
            'scaling', 'linear', 'controlSteps', i);
    end
    
    % ---- Objective functions using observed states ----
    obj_func_bfgs = @(m1,s1,sch1,m2,s2,sch2,sc,obj_prior, deriv,t,s1_in,s2_in) ...
        matchObservedTwinModelsOW(m1,s1,sch1,m2,s2,sch2,sc,obj_prior,...
        'computePartials',deriv,'tstep',t,'state1',s1_in,...
        'state2',s2_in,'from_states',false,weighting{:});
    
    obj_bfgs = @(p) evaluateMatchTwinModels(p, obj_func_bfgs, setup_BhpRate, ...
                                             setup_RateBhp, params_BhpRate, ...
                                             statesObserved, prior);
    
    obj_func_lm = @(m1,s1,sch1,m2,s2,sch2,sc,obj_prior, deriv,t,s1_in,s2_in) ...
        matchObservedTwinModelsOW(m1,s1,sch1,m2,s2,sch2,sc, obj_prior,...
        'computePartials',deriv,'tstep',t,'state1',s1_in,...
        'state2',s2_in,'from_states',false,weighting{:}, 'mismatchSum', false);
    
    obj_lm = @(p) evaluateMatchTwinModelsSummands(p, obj_func_lm, setup_BhpRate, ...
                                                  setup_RateBhp, params_RateBhp, ...
                                                  statesObserved, prior);
    
    % ---- BFGS Stage ----
    use_bfgs = true;
    bfgs_iterations = 5;
    p0 = getScaledParameterVector(setup_BhpRate, params_BhpRate);
    if use_bfgs
        fprintf('BFGS optimization...\n');
        [~, p_opt, ~] = unitBoxBFGS(p0, obj_bfgs, ...
            'maxIt', bfgs_iterations, 'objChangeTol', 1e-7, 'gradTol', 1e-7);
    else
        p_opt = p0;
    end
    
    % ---- Levenberg-Marquardt with history recording ----
    fprintf('Levenberg-Marquardt optimization...\n');
    [~, p_final, info] = unitBoxLM(p_opt, obj_lm, 'maxIt', 30);
    obj_history = info.val;   % scalar objective per iteration
    p_opt_final = p_final;
    
    % Store history
    objHistories{iNoise} = [(1:length(obj_history)), obj_history];
    
    % Update final model from optimized parameters
    setup_final = updateSetupFromScaledParameters(setup_BhpRate, params_BhpRate, p_opt_final);
    finalModels{iNoise} = setup_final.model;
    finalParams{iNoise} = p_opt_final;
    
    % Simulate final optimized model to get well solutions
    [wellSols_opt, ~] = simulateScheduleAD(state0, setup_final.model, setup_final.schedule);
    
    % Extract recovered BHP from optimized schedule
    recBhpP1 = zeros(numel(setup_final.schedule.control),1);
    recBhpP2 = zeros(numel(setup_final.schedule.control),1);
    for iCtrl = 1:numel(setup_final.schedule.control)
        recBhpP1(iCtrl) = setup_final.schedule.control(iCtrl).W(3).val / barsa;
        recBhpP2(iCtrl) = setup_final.schedule.control(iCtrl).W(4).val / barsa;
    end
    recoveredBhpP1{iNoise} = recBhpP1;
    recoveredBhpP2{iNoise} = recBhpP2;
    
    % Extract recovered water and oil rates
    recWaterP1 = zeros(numel(wellSols_opt),1);
    recWaterP2 = zeros(numel(wellSols_opt),1);
    recOilP1   = zeros(numel(wellSols_opt),1);
    recOilP2   = zeros(numel(wellSols_opt),1);
    for t = 1:numel(wellSols_opt)
        recWaterP1(t) = abs(wellSols_opt{t}(3).qWs);
        recWaterP2(t) = abs(wellSols_opt{t}(4).qWs);
        recOilP1(t)   = abs(wellSols_opt{t}(3).qOs);
        recOilP2(t)   = abs(wellSols_opt{t}(4).qOs);
    end
    recoveredWaterRateP1{iNoise} = recWaterP1;
    recoveredWaterRateP2{iNoise} = recWaterP2;
    recoveredOilRateP1{iNoise}   = recOilP1;
    recoveredOilRateP2{iNoise}   = recOilP2;
    
    fprintf('Finished noise level %.1f%%, final objective = %.6e\n', relNoise*100, obj_history(end));
end

%% (7) Final Figures – All Noise Levels

time = cumsum([0; scheduleRef.step.val(1:end-1)]) / day;
colors = lines(numNoise);

% ----- Figure 4: BHP for Well P1 and Well P2 (2×1) with all noise levels -----
figure(4);
clf;
subplot(2,1,1);
plot(time, trueBhpP1, 'k-', 'LineWidth', 2, 'DisplayName', 'True BHP'); hold on;
for i = 1:numNoise
    plot(time, recoveredBhpP1{i}, '--', 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('Noise %.1f%%', noiseLevels(i)*100));
end
xlabel('Time (days)'); ylabel('BHP (bar)');
legend('Location','eastoutside'); title('Well P1 - BHP reconstruction'); grid on;

subplot(2,1,2);
plot(time, trueBhpP2, 'k-', 'LineWidth', 2, 'DisplayName', 'True BHP'); hold on;
for i = 1:numNoise
    plot(time, recoveredBhpP2{i}, '--', 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('Noise %.1f%%', noiseLevels(i)*100));
end
xlabel('Time (days)'); ylabel('BHP (bar)');
legend('Location','eastoutside'); title('Well P2 - BHP reconstruction'); grid on;

% ----- Figure 5: Water and oil rates (2×2) with all noise levels -----
figure(5);
clf;
% Well P1 water rate
subplot(2,2,1);
plot(time, trueWaterRateP1, 'k-', 'LineWidth', 2, 'DisplayName', 'True water rate'); hold on;
for i = 1:numNoise
    plot(time, recoveredWaterRateP1{i}, '--', 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('Noise %.1f%%', noiseLevels(i)*100));
end
xlabel('Time (days)'); ylabel('Water rate (m^3/day)');
legend('Location','eastoutside'); title('Well P1 - Water rate'); grid on;

% Well P1 oil rate
subplot(2,2,2);
plot(time, trueOilRateP1, 'k-', 'LineWidth', 2, 'DisplayName', 'True oil rate'); hold on;
for i = 1:numNoise
    plot(time, recoveredOilRateP1{i}, '--', 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('Noise %.1f%%', noiseLevels(i)*100));
end
xlabel('Time (days)'); ylabel('Oil rate (m^3/day)');
legend('Location','eastoutside'); title('Well P1 - Oil rate'); grid on;

% Well P2 water rate
subplot(2,2,3);
plot(time, trueWaterRateP2, 'k-', 'LineWidth', 2, 'DisplayName', 'True water rate'); hold on;
for i = 1:numNoise
    plot(time, recoveredWaterRateP2{i}, '--', 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('Noise %.1f%%', noiseLevels(i)*100));
end
xlabel('Time (days)'); ylabel('Water rate (m^3/day)');
legend('Location','eastoutside'); title('Well P2 - Water rate'); grid on;

% Well P2 oil rate
subplot(2,2,4);
plot(time, trueOilRateP2, 'k-', 'LineWidth', 2, 'DisplayName', 'True oil rate'); hold on;
for i = 1:numNoise
    plot(time, recoveredOilRateP2{i}, '--', 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('Noise %.1f%%', noiseLevels(i)*100));
end
xlabel('Time (days)'); ylabel('Oil rate (m^3/day)');
legend('Location','eastoutside'); title('Well P2 - Oil rate'); grid on;

% ----- Figure 6: RMSE convergence vs noise level (2×2) -----
figure(6);
clf;
% Compute RMSE for each noise level
rmse_bhp_p1 = zeros(numNoise,1);
rmse_bhp_p2 = zeros(numNoise,1);
rmse_water_p1 = zeros(numNoise,1);
rmse_water_p2 = zeros(numNoise,1);
rmse_oil_p1   = zeros(numNoise,1);
rmse_oil_p2   = zeros(numNoise,1);
for i = 1:numNoise
    rmse_bhp_p1(i) = sqrt(mean((recoveredBhpP1{i} - trueBhpP1).^2));
    rmse_bhp_p2(i) = sqrt(mean((recoveredBhpP2{i} - trueBhpP2).^2));
    rmse_water_p1(i) = sqrt(mean((recoveredWaterRateP1{i} - trueWaterRateP1).^2));
    rmse_water_p2(i) = sqrt(mean((recoveredWaterRateP2{i} - trueWaterRateP2).^2));
    rmse_oil_p1(i)   = sqrt(mean((recoveredOilRateP1{i} - trueOilRateP1).^2));
    rmse_oil_p2(i)   = sqrt(mean((recoveredOilRateP2{i} - trueOilRateP2).^2));
end

% Subplot 1: BHP RMSE
subplot(2,2,1);
bar(noiseLevels*100, [rmse_bhp_p1, rmse_bhp_p2]);
xlabel('Noise level (%)'); ylabel('RMSE (bar)');
legend('Well P1', 'Well P2', 'Location','best');
title('BHP reconstruction error'); grid on;

% Subplot 2: Water rate RMSE
subplot(2,2,2);
bar(noiseLevels*100, [rmse_water_p1, rmse_water_p2]);
xlabel('Noise level (%)'); ylabel('RMSE (m^3/day)');
legend('Well P1', 'Well P2', 'Location','best');
title('Water rate reconstruction error'); grid on;

% Subplot 3: Oil rate RMSE
subplot(2,2,3);
bar(noiseLevels*100, [rmse_oil_p1, rmse_oil_p2]);
xlabel('Noise level (%)'); ylabel('RMSE (m^3/day)');
legend('Well P1', 'Well P2', 'Location','best');
title('Oil rate reconstruction error'); grid on;

% Subplot 4: Objective convergence (all noise levels)
subplot(2,2,4);
hold on;
for i = 1:numNoise
    hist = objHistories{i};
    semilogy(hist(:,1), hist(:,2), 'Color', colors(i,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Noise %.1f%%', noiseLevels(i)*100));
end
xlabel('Iteration'); ylabel('Objective function');
legend('Location','best'); title('Convergence of objective');
grid on;

fprintf('\nAll figures generated.\n');
fprintf('Figure 4: BHP reconstruction for all noise levels.\n');
fprintf('Figure 5: Water & oil rates reconstruction for all noise levels.\n');
fprintf('Figure 6: RMSE vs noise level and objective convergence.\n');

fprintf('\nAll tests completed. Figures show reconstruction of production well BHP and rates.\n');