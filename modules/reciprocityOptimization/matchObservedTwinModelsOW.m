function obj = matchObservedTwinModelsOW(model_1, states_1, schedule_1, model_2, states_2, schedule_2, observed, prior_dist, varargin)
% Compute mismatch between twin reservoir models with improved organization and scaling
%
% Mismatch types:
%   'leastSquares' - J = sum(Δt/2 * [ωq^{-1}*(δ(ρq))^2 + ωp*(δp)^2])
%   'energyWell'   - J = sum(Δt/2 * [(ρq)_A - (ρq)_B] * [p_A - p_B])
%   'reciprocity'  - J = sum(Δt/2 * [(ρq)_A * p_B - (ρq)_B * p_A])
%
% Regularization terms (for leastSquares):
%   'pressure'     - K * ρλ * (δp)^2  (semi-norm weighted by properties)
%   'massAccum'    - |∂t δ(ρ s)|^2
%   'saturation'   - |δ s|^2

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
%}

% Enhanced options with better default scaling
opt = struct('WaterRateWeight',     [], ...
    'OilRateWeight',       [], ...
    'BHPWeight',           [], ...
    'ComputePartials',     false, ...
    'tStep',               [], ...
    'state1',              [], ...
    'state2',              [], ...
    'from_states',         false, ...
    'matchOnlyProducers',  false, ...
    'mismatchSum',         true, ...
    'accumulateWells',     [], ...
    'accumulateTypes',     [], ...
    'regBeta',             1e-3, ...       % Regularization weight β (smaller default)
    'regPressure',         true, ...        % K * ρλ * (δp)^2 term
    'regMassAccum',        false, ...       % |∂t δ(ρ s)|^2 term
    'regSaturation',       false, ...       % |δ s|^2 term
    'scalingMethod',       'dynamic', ...
    'rho_ref_water',       1000, ...        % Reference density for water (kg/m³)
    'rho_ref_oil',         850);            % Reference density for oil (kg/m³)

opt = merge_options(opt, varargin{:});

% Time step handling
dts = schedule_1.step.val;
totTime = sum(dts);

if isempty(opt.tStep)
    tSteps = (1:numel(dts))';
else
    tSteps = opt.tStep;
    dts = dts(opt.tStep);
end

obj = repmat({[]}, numel(tSteps), 1);

for step = 1:numel(tSteps)
    current_step = tSteps(step);
    dt = dts(step);

    % Get all data in a more organized way
    [modelData1, modelData2, obsData] = getAllModelData(...
        model_1, states_1, ...
        model_2, states_2, ...
        observed, current_step, opt);

    % Calculate weights based on observed data scaling with reference densities
    [ww, wo, wp] = getObservedWeights(obsData.qWs, obsData.qOs, obsData.bhp, ...
        opt.rho_ref_water, opt.rho_ref_oil, opt);

    % Select mismatch type
    mismatchType = 'leastSquares';  % Default
    
    % Compute mismatch components based on selected type
    mismatchComponents = computeMismatchComponents(...
        modelData1, modelData2, schedule_1, obsData, ...
        ww, wo, wp, dt, totTime, mismatchType, opt);

    % Add regularization terms for least squares
    if strcmp(mismatchType, 'leastSquares') && opt.regBeta > 0
        regComponents = computeRegularizationComponents(...
            model_1, modelData1, modelData2, dt, totTime, opt);
        well_cells = vertcat(schedule_1.control(1).W.cells);
        regOnlywells = true;
        if regOnlywells
            for i = length(regComponents)
                regComponents_wells{i} = opt.regBeta.*regComponents{i}(well_cells);
            end
        end
        aa{1} = opt.regBeta.*sum(regComponents{1});
        mismatchComponents = [mismatchComponents; regComponents_wells;aa];
        
%         % Add regularization components to mismatch
%         for i = 1:numel(regComponents)
%             mismatchComponents{end+1} = regComponents{i};
%         end
    end

    % Handle accumulation options
    if ~isempty(opt.accumulateTypes)
        pt = opt.accumulateTypes;
        tmp = num2cell(zeros(1, max(pt)));
        for k = 1:numel(mismatchComponents)
            if pt(k) > 0
                tmp{pt(k)} = tmp{pt(k)} + mismatchComponents{k};
            end
        end
        mismatchComponents = tmp;
    end

    if ~isempty(opt.accumulateWells)
        pw = opt.accumulateWells;
        M = sparse(pw(pw>0), find(pw), 1);
        mismatchComponents = cellfun(@(x) M*x, mismatchComponents, 'UniformOutput', false);
    end

    % Return either summed or component form
    if opt.mismatchSum
        obj{step} = sum(vertcat(mismatchComponents{:}));
    else
        obj{step} = vertcat(mismatchComponents{:});
    end

    % Add prior term if provided
    if ~isempty(prior_dist) && opt.mismatchSum
        obj{step} = obj{step} + 0 * sum(prior_dist{step}.^2);
    end
end

% If single step, return directly instead of cell array
if numel(obj) == 1
    obj = obj{1};
end
end

%% Helper Functions

function [model1, model2, obs] = getAllModelData(model_1, states_1, ...
    model_2, states_2, ...
    observed, current_step, opt)
% Organized data extraction for both models and observed data

% Observed data
sol_obs = observed{current_step};
nw = numel(sol_obs.wellSol);
obs.status = vertcat(sol_obs.wellSol.status);
obs.qWs = vertcatIfPresent(sol_obs.wellSol, 'qWs', nw);
obs.qOs = vertcatIfPresent(sol_obs.wellSol, 'qOs', nw);
obs.bhp = vertcatIfPresent(sol_obs.wellSol, 'bhp', nw);
obs.sign = vertcat(sol_obs.wellSol.sign);

% Model 1 data
if opt.ComputePartials && opt.from_states
    state_1 = model_1.getStateAD(states_1{current_step}, true);
elseif opt.ComputePartials
    state_1 = opt.state1;
else
    state_1 = states_1{current_step};
end

model1 = extractModelData(model_1, state_1, opt.ComputePartials);

% Model 2 data
if opt.ComputePartials && opt.from_states
    state_2 = model_2.getStateAD(states_2{current_step}, true);
elseif opt.ComputePartials
    state_2 = opt.state2;
else
    state_2 = states_2{current_step};
end

model2 = extractModelData(model_2, state_2, opt.ComputePartials);

% Handle inactive wells consistently
[model1, model2, obs] = handleInactiveWells(model1, model2, obs);
end

function modelData = extractModelData(model, state, isAD)
% Extract all relevant data from a model state
nw = numel(state.wellSol);

modelData.status = vertcat(state.wellSol.status);

if isAD && isa(state.pressure, 'ADI')
    modelData.pressure = model.getProp(state, 'pressure');
    modelData.saturation = model.getProp(state, 's');
    modelData.qWs = model.FacilityModel.getProp(state, 'qWs');
    modelData.qOs = model.FacilityModel.getProp(state, 'qOs');
    modelData.bhp = model.FacilityModel.getProp(state, 'bhp');
    modelData.wellSol = state.wellSol;
    modelData.totalFlux = state.FluxDisc.ComponentTotalFlux;
    modelData.totalMass = state.FlowProps.ComponentTotalMass;
    modelData.Mobility = state.FlowProps.Mobility;
    modelData.FaceComponentMobility = state.FluxDisc.FaceComponentMobility;
    modelData.density = state.PVTProps.Density;
else
    modelData.pressure = state.pressure;
    modelData.saturation = state.s;
    modelData.totalFlux = state.FluxDisc.ComponentTotalFlux;
    modelData.totalMass = state.FlowProps.ComponentTotalMass;
    modelData.Mobility = state.FlowProps.Mobility;
    modelData.FaceComponentMobility = state.FluxDisc.FaceComponentMobility;
    modelData.density = state.PVTProps.Density;
    modelData.qWs = vertcatIfPresent(state.wellSol, 'qWs', nw);
    modelData.qOs = vertcatIfPresent(state.wellSol, 'qOs', nw);
    modelData.bhp = vertcatIfPresent(state.wellSol, 'bhp', nw);
    modelData.wellSol = state.wellSol;
end

% Ensure pressure is cell array for consistency
if ~iscell(modelData.pressure)
    modelData.pressure = {modelData.pressure};
    modelData.saturation = {modelData.saturation};
end
end

function [model1, model2, obs] = handleInactiveWells(model1, model2, obs)
% Handle inactive wells consistently across models and observations
if ~all(model1.status) || ~all(model2.status) || ~all(obs.status)
    [model1.bhp, obs.bhp] = expandToFull(model1.bhp, obs.bhp, model1.status, obs.status, true);
    [model1.qWs, obs.qWs] = expandToFull(model1.qWs, obs.qWs, model1.status, obs.status, false);
    [model1.qOs, obs.qOs] = expandToFull(model1.qOs, obs.qOs, model1.status, obs.status, false);

    [model2.bhp, obs.bhp] = expandToFull(model2.bhp, obs.bhp, model2.status, obs.status, true);
    [model2.qWs, obs.qWs] = expandToFull(model2.qWs, obs.qWs, model2.status, obs.status, false);
    [model2.qOs, obs.qOs] = expandToFull(model2.qOs, obs.qOs, model2.status, obs.status, false);
end
end

function [ww, wo, wp] = getObservedWeights(qWs, qOs, bhp, rho_ref_w, rho_ref_o, opt)
% Calculate weights based on observed data scaling as in equation:
% ωq = ωp = ||q_obs||_ℓ2 / ||p_obs||_ℓ2
% But now properly scaled for mass rate mismatch: (ρq) vs q_obs

% If weights are provided explicitly, use them
if ~isempty(opt.WaterRateWeight) && ~isempty(opt.OilRateWeight) && ~isempty(opt.BHPWeight)
    ww = opt.WaterRateWeight;
    wo = opt.OilRateWeight;
    wp = opt.BHPWeight;
    return;
end

% Compute ℓ2 norms of observed volumetric rates
norm_qW = sqrt(sum(qWs.^2));
norm_qO = sqrt(sum(qOs.^2));
norm_p = sqrt(sum(bhp.^2));

% Avoid division by zero
if norm_p == 0
    scale = 1;
else
    scale = 1/norm_p;
end

% Set weights according to corrected scaling:
% For mass rate mismatch (ρq), the observed data is q (volumetric)
% So ωq should be: ρ_ref * ||q_obs|| / ||p_obs||
if isempty(opt.WaterRateWeight)
    if norm_qW == 0
        ww = 0;
    else
        ww = rho_ref_w * norm_qW * scale;  % Now includes ρ for water
    end
end

if isempty(opt.OilRateWeight)
    if norm_qO == 0
        wo = 0;
    else
        wo = rho_ref_o * norm_qO * scale;  % Now includes ρ for oil
    end
end

if isempty(opt.BHPWeight)
    if norm_p == 0
        wp = 0;
    else
        wp = norm_p * scale;  % Unchanged (pressure mismatch)
    end
end

% Ensure weights are reasonable (avoid extreme values)
ww = min(max(ww, 1e-6), 1e6);
wo = min(max(wo, 1e-6), 1e6);
wp = min(max(wp, 1e-6), 1e6);
end

function v = vertcatIfPresent(sol, fn, nw)
if isfield(sol, fn)
    v = vertcat(sol.(fn));
    assert(numel(v) == nw);
    v = v(vertcat(sol.status));
else
    v = zeros(nnz(sol.status), 1);
end
end

function [v, v_obs] = expandToFull(v, v_obs, status, status_obs, setToZero)
tmp = zeros(size(status));
if isa(v, 'ADI')
    tmp = double2ADI(tmp, v);
end
tmp(status) = v;
v = tmp;

tmp = zeros(size(status));
tmp(status_obs) = v_obs;
v_obs = tmp;

if setToZero
    ix = status ~= status_obs;
    v(ix) = 0;
    v_obs(ix) = 0;
end
end

%% Main Mismatch Component Function

function components = computeMismatchComponents(model1, model2, schedule, obs, ww, wo, wp, dt, totTime, mismatchType, opt)
% Compute mismatch components for all wells and phases based on selected type

matchCases = true(numel(obs.qWs), 1);
if opt.matchOnlyProducers
    matchCases = (obs.sign < 0);
end

% Get well cells
well_cells = vertcat(schedule.control(1).W.cells);

if isempty(well_cells)
    components = {0};
    return;
end

% Get densities at wells - handle each model independently
if iscell(model1.density)
    rho1 = model1.density{1};
else
    rho1 = model1.density(:,1);
end

if iscell(model2.density)
    rho2 = model2.density{1};
else
    rho2 = model2.density(:,1);
end

% Get pressures - handle each model independently
if iscell(model1.pressure)
    p1_all = model1.pressure{1};
else
    p1_all = model1.pressure;
end

if iscell(model2.pressure)
    p2_all = model2.pressure{1};
else
    p2_all = model2.pressure;
end

% Extract pressures at wells
p1_well = p1_all(well_cells);
p2_well = p2_all(well_cells);
rho1_well = rho1(well_cells);
rho2_well = rho2(well_cells);

% Get well rates
qWs1 = model1.qWs;
qOs1 = model1.qOs;
qWs2 = model2.qWs;
qOs2 = model2.qOs;

% Apply matching cases filter
qWs1 = qWs1 .* matchCases;
qOs1 = qOs1 .* matchCases;
qWs2 = qWs2 .* matchCases;
qOs2 = qOs2 .* matchCases;

% Compute mass rates (ρq) for each phase and well
massRate1_w = rho1_well .* qWs1;
massRate1_o = rho1_well .* qOs1;
massRate2_w = rho2_well .* qWs2;
massRate2_o = rho2_well .* qOs2;

fac = dt / (2 * totTime);
nw = numel(well_cells);

% Ensure weights are per-well vectors (if scalar, expand to vector)
if isscalar(ww)
    ww = ww * ones(nw, 1);
end
if isscalar(wo)
    wo = wo * ones(nw, 1);
end
if isscalar(wp)
    wp = wp * ones(nw, 1);
end

switch mismatchType
    case 'leastSquares'
        % J = sum(Δt/2 * [ωq^{-1}*(δ(ρq))^2 + ωp*(δp)^2])
        components = cell(nw * 3, 1); % water, oil, pressure for each well
        
        for w = 1:nw
            % Water rate component: ωq^{-1} * (δ(ρq))^2
            delta_massRate_w = massRate1_w(w) - massRate2_w(w);
            if ww(w) > 0
                components{w} = fac * (delta_massRate_w^2) / ww(w);
            else
                components{w} = 0;
            end
            
            % Oil rate component: ωq^{-1} * (δ(ρq))^2
            delta_massRate_o = massRate1_o(w) - massRate2_o(w);
            if wo(w) > 0
                components{nw + w} = fac * (delta_massRate_o^2) / wo(w);
            else
                components{nw + w} = 0;
            end
            
            % Pressure component: ωp * (δp)^2
            delta_p = p1_well(w) - p2_well(w);
            components{2*nw + w} = fac * wp(w) * (delta_p^2);
        end
        
    case 'energyWell'
        % J = sum(Δt/2 * [(ρq)_A - (ρq)_B] * [p_A - p_B])
        components = cell(nw * 2, 1); % water and oil for each well
        
        for w = 1:nw
            delta_p = p1_well(w) - p2_well(w);
            
            % Water component
            delta_massRate_w = massRate1_w(w) - massRate2_w(w);
            components{w} = abs(fac * delta_massRate_w * delta_p);
            
            % Oil component
            delta_massRate_o = massRate1_o(w) - massRate2_o(w);
            components{nw + w} = abs(fac * delta_massRate_o * delta_p);
        end
        
    case 'reciprocity'
        % J = sum(Δt/2 * [(ρq)_A * p_B - (ρq)_B * p_A])
        components = cell(nw * 2, 1); % water and oil for each well
        
        for w = 1:nw
            % Water component
            components{w} = fac * abs(massRate1_w(w) * p2_well(w) - massRate2_w(w) * p1_well(w));
            
            % Oil component
            components{nw + w} = fac * abs(massRate1_o(w) * p2_well(w) - massRate2_o(w) * p1_well(w));
        end
        
    otherwise
        error('Unknown mismatch type: %s', mismatchType);
end
end

function regComponents = computeRegularizationComponents(model, model1, model2, dt, totTime, opt)
% Compute regularization terms from the paper:
%   regPressure  - K * ρλ * (δp)^2  (semi-norm weighted by properties)
%   regMassAccum - |∂t δ(ρ s)|^2
%   regSaturation - |δ s|^2

regComponents = {};
idx = 1;
fac = dt / (2 * totTime);  % Consistent with mismatch scaling

% Get base pressure fields - handle each model independently
if iscell(model1.pressure)
    p1 = model1.pressure{1};
else
    p1 = model1.pressure;
end

if iscell(model2.pressure)
    p2 = model2.pressure{1};
else
    p2 = model2.pressure;
end

delta_p = p1 - p2;

% Get number of phases from model1 (assuming it's the reference)
if iscell(model1.saturation)
    nph = numel(model1.saturation);
else
    nph = size(model1.saturation, 2);
end

% Pressure regularization: K * ρλ * (δp)^2 (semi-norm)
if opt.regPressure
    K = model.rock.perm(:,1);  % Permeability (assuming isotropic)
    regPressure = zeros(size(p1));
    
    for ph = 1:nph
        % Get density for this phase from model1
        if iscell(model1.density)
            rho1 = model1.density{ph};
        else
            rho1 = model1.density(:,ph);
        end
        
        % Get density for this phase from model2
        if iscell(model2.density)
            if ph <= numel(model2.density)
                rho2 = model2.density{ph};
            else
                rho2 = rho1;  % Fallback to model1 if model2 missing
            end
        else
            if ph <= size(model2.density, 2)
                rho2 = model2.density(:,ph);
            else
                rho2 = rho1;  % Fallback to model1 if model2 missing
            end
        end
        rho_avg = (rho1 + rho2) / 2;  % Average density
        
        % Get mobility for this phase from model1
        if iscell(model1.Mobility)
            mob1 = model1.Mobility{ph};
        else
            mob1 = model1.Mobility(:,ph);
        end
        
        % Get mobility for this phase from model2
        if iscell(model2.Mobility)
            if ph <= numel(model2.Mobility)
                mob2 = model2.Mobility{ph};
            else
                mob2 = mob1;  % Fallback to model1 if model2 missing
            end
        else
            if ph <= size(model2.Mobility, 2)
                mob2 = model2.Mobility(:,ph);
            else
                mob2 = mob1;  % Fallback to model1 if model2 missing
            end
        end
        mob_avg = (mob1 + mob2) / 2;  % Average mobility
        
        % Semi-norm: ρλ (p1-p2)^2
        regPressure = regPressure +  mob_avg .* delta_p.^2;
    end
    
    % Scale by permeability and add to components
    regComponents{idx} = fac * opt.regBeta .* K .* regPressure;
    idx = idx + 1;
end

% Mass accumulation regularization: |∂t δ(ρ s)|^2
if opt.regMassAccum
    regMassAccum = zeros(size(p1));
    
    for ph = 1:nph
        % Get density for this phase from model1
        if iscell(model1.density)
            rho1 = model1.density{ph};
        else
            rho1 = model1.density(:,ph);
        end
        
        % Get density for this phase from model2
        if iscell(model2.density)
            if ph <= numel(model2.density)
                rho2 = model2.density{ph};
            else
                rho2 = rho1;  % Fallback to model1 if model2 missing
            end
        else
            if ph <= size(model2.density, 2)
                rho2 = model2.density(:,ph);
            else
                rho2 = rho1;  % Fallback to model1 if model2 missing
            end
        end
        
        % Get saturation for this phase from model1
        if iscell(model1.saturation)
            s1 = model1.saturation{ph};
        else
            s1 = model1.saturation(:,ph);
        end
        
        % Get saturation for this phase from model2
        if iscell(model2.saturation)
            if ph <= numel(model2.saturation)
                s2 = model2.saturation{ph};
            else
                s2 = s1;  % Fallback to model1 if model2 missing
            end
        else
            if ph <= size(model2.saturation, 2)
                s2 = model2.saturation(:,ph);
            else
                s2 = s1;  % Fallback to model1 if model2 missing
            end
        end
        
        % Mass difference between simulators
        delta_mass = rho1 .* s1 - rho2 .* s2;
        
        % Time derivative approximation
        regMassAccum = regMassAccum + (delta_mass / dt).^2;
    end
    
    regComponents{idx} = fac * opt.regBeta .* regMassAccum;
    idx = idx + 1;
end

% Saturation regularization: |δ s|^2
if opt.regSaturation
    regSaturation = zeros(size(p1));
    
    for ph = 1:nph
        % Get saturation for this phase from model1
        if iscell(model1.saturation)
            s1 = model1.saturation{ph};
        else
            s1 = model1.saturation(:,ph);
        end
        
        % Get saturation for this phase from model2
        if iscell(model2.saturation)
            if ph <= numel(model2.saturation)
                s2 = model2.saturation{ph};
            else
                s2 = s1;  % Fallback to model1 if model2 missing
            end
        else
            if ph <= size(model2.saturation, 2)
                s2 = model2.saturation(:,ph);
            else
                s2 = s1;  % Fallback to model1 if model2 missing
            end
        end
        
        delta_s = s1 - s2;
        regSaturation = regSaturation + delta_s.^2;
    end
    
    regComponents{idx} = fac * opt.regBeta .* regSaturation;
    idx = idx + 1;
end

% If no regularization terms selected, return empty
if isempty(regComponents)
    regComponents = {0};
end
end