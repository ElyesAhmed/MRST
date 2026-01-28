classdef ComponentPhaseMolecularDiffFlux < StateFunction
    % ComponentPhaseMolecularDiffFlux - Molecular diffusion flux for each component
    %
    % Uses Quirk-Millington tortuosity model with Blanc's law for gas phase
    % binary diffusion and literature values for liquid phase.

    methods
        function gp = ComponentPhaseMolecularDiffFlux(model, varargin)
            gp@StateFunction(model);
            gp = setupDependencies(gp, model);
            gp.label = 'J_{i,\alpha}';
        end

        function J = evaluateOnDomain(prop, model, state)
            % Early return if no compositional data
            if ~isfield(state, 'x')
                ncomp = model.getNumberOfComponents;
                nph = model.getNumberOfPhases;
                J = cell(ncomp, nph);
                J = cellfun(@(x) 0, J, 'UniformOutput', false);
                return;
            end
            J = calculateMolecularDiffusionFluxes(prop, model, state);
        end
    end
end

%% Setup and Calculation Functions
function gp = setupDependencies(gp, model)
gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
gp = gp.dependsOn('PoreVolume', 'PVTPropertyFunctions');
gp = gp.dependsOn('s', 'state');
gp = gp.dependsOn('x', 'state');
gp = gp.dependsOn('y', 'state');
gp = gp.dependsOn('pressure', 'state');
gp = gp.dependsOn('T', 'state');
end

function J = calculateMolecularDiffusionFluxes(prop, model, state)
ncomp = model.getNumberOfComponents;
nph = model.getNumberOfPhases;
nm = model.getPhaseNames();

% Initialize
J = cell(ncomp, nph);
[J{:}] = deal(0);

% Get diffusion parameters
[mol_diff, param_LJ, Dij] = getDiffusionParameters(model);

% Get required properties
rho = prop.getEvaluatedExternals(model, state, 'Density');
pvtprops = model.PVTPropertyFunctions;
pv = pvtprops.get(prop, state, 'PoreVolume');
[p, T] = model.getProps(state, 'pressure', 'temperature');
[xc, yc] = getMoleFractions(model, state);
% Gas diffusion coefficient base
p_atm = convertTo(value(p), atm);
coeff1 = T.^(1.5) ./ p_atm;

% Phase indices
L_ix = model.getLiquidIndex();
V_ix = model.getVaporIndex();

% Calculate fluxes
for c = 1:ncomp

    for ph = 1:nph
        s = model.getProp(state, ['s', nm(ph)]);
        if ph == L_ix
            J{c, ph} = calculateLiquidDiffusionFlux(...
                model, s, rho{ph}, xc{c}, mol_diff(c, ph), pv);

        elseif ph == V_ix
            J{c, ph} = calculateGasDiffusionFlux(...
                model, state, s, rho{ph}, yc{c}, coeff1, Dij(c, :), pv, c);
        end
    end
end
end

%% Diffusion Coefficient Functions
function [mol_diff, param_LJ, Dij] = getDiffusionParameters(model)
% Get component names
namecp = model.compFluid.names();

% Component database
componentNames = {'H2', 'C1', 'CO2', 'H2O', 'N2', 'C2', 'C3', 'H2S', 'NC4', ...
    'water', 'oil', 'gas', 'methane', 'ethane', 'propane', 'butane'};

indices = createComponentIndices(namecp, componentNames);
[mol_diff, param_LJ] = loadDiffusionCoefficients(indices, componentNames);
Dij = calculateBinaryDiffusionMatrix(model, param_LJ);
end

function indices = createComponentIndices(namecp, componentNames)
indices = struct();
for i = 1:numel(componentNames)
    name = componentNames{i};
    idx = find(strcmp(namecp, name));
    if ~isempty(idx)
        indices.(lower(name)) = idx(1);
    end
end
end

function [mol_diff, param_LJ] = loadDiffusionCoefficients(indices, componentNames)
% Default values
DEFAULT_LIQUID = 1e-9;
DEFAULT_GAS = 1e-5;
DEFAULT_LJ = [3.5, 150.0];

% Determine number of components
if isempty(fieldnames(indices))
    error('No component names matched. Check component naming.');
end
ncomp = max(cell2mat(struct2cell(indices)));

mol_diff = DEFAULT_LIQUID * ones(ncomp, 2);
param_LJ = repmat(DEFAULT_LJ, ncomp, 1);

% Database (40°C) [liquid, gas] in m²/s
coeffs = struct(...
    'h2',       [6.44e-9, 6.1e-5], ...
    'c1',       [2.15e-9, 1.6e-5], ...
    'methane',  [2.15e-9, 1.6e-5], ...
    'co2',      [2.72e-9, 1.4e-5], ...
    'h2o',      [3.29e-9, 1.5e-5], ...
    'water',    [3.29e-9, 1.5e-5], ...
    'n2',       [2.86e-9, 1.8e-5], ...
    'c2',       [1.72e-9, 2.5e-5], ...
    'ethane',   [1.72e-9, 2.5e-5], ...
    'c3',       [1.43e-9, 2.2e-5], ...
    'propane',  [1.43e-9, 2.2e-5], ...
    'h2s',      [2.15e-9, 2.2e-5], ...
    'nc4',      [1.15e-9, 1.9e-5], ...
    'butane',   [1.15e-9, 1.9e-5], ...
    'oil',      [1.0e-9,  1.5e-5], ...
    'gas',      [0.5e-9,  2.0e-5]);

coeffs_LJ = struct(...
    'h2',       [2.92,   59.7], ...
    'c1',       [3.758, 148.6], ...
    'methane',  [3.758, 148.6], ...
    'co2',      [3.996, 195.2], ...
    'h2o',      [2.641, 809.1], ...
    'water',    [2.641, 809.1], ...
    'n2',       [3.798,  71.4], ...
    'c2',       [4.443, 215.7], ...
    'ethane',   [4.443, 215.7], ...
    'c3',       [5.118, 237.1], ...
    'propane',  [5.118, 237.1], ...
    'h2s',      [3.60,  301.0], ...
    'nc4',      [5.206, 289.5], ...
    'butane',   [5.206, 289.5], ...
    'oil',      [5.0,   400.0], ...
    'gas',      [3.5,   100.0]);

% Assign known coefficients
fields = fieldnames(indices);
for i = 1:numel(fields)
    comp = fields{i};
    idx = indices.(comp);
    if isfield(coeffs, comp)
        mol_diff(idx, :) = coeffs.(comp);
        param_LJ(idx, :) = coeffs_LJ.(comp);
    end
end
end

function Dij = calculateBinaryDiffusionMatrix(model, param_LJ)
ncomp = size(param_LJ, 1);
Molmass = 1.e3 .* model.compFluid.molarMass;
SigLJ = param_LJ(:, 1);

Dij = zeros(ncomp);
for i = 1:ncomp
    for j = 1:ncomp
        sqrtMij = sqrt(2 * Molmass(i) * Molmass(j) / (Molmass(i) + Molmass(j)));
        Sigij2 = 0.25 * (SigLJ(i) + SigLJ(j))^2;
        Dij(i, j) = 1.e-4 * 0.001858 / (sqrtMij * Sigij2); % m²/s
    end
end
end

%% Flux Calculation Functions
function [xc, yc] = getMoleFractions(model, state)
% For black-oil models
if isfield(state, 'rs') || isfield(state, 'rv')
    error('Black-oil diffusion with rs/rv is under development');
end

[x, y] = model.getProps(state, 'x', 'y');
% Compositional models
if iscell(x)
    xc = x;
else        
    xc = mat2cell(x, size(x, 1), ones(1, size(x, 2)));
end

if iscell(y)
    yc = y;
else
    yc = mat2cell(y, size(y, 1), ones(1, size(y, 2)));
end
end

function flux = calculateLiquidDiffusionFlux(model, s, rho, xc, D_mol, poro)
op = model.operators;

% Millington-Quirk tortuosity
tau_mq = (s .* poro).^(7/3) .* poro.^(-2);

% Effective diffusivity
D_eff = op.faceAvg(s .* rho .* D_mol .* tau_mq .* poro);

% Diffusion flux
flux = -D_eff .* op.Grad(xc);
end

function flux = calculateGasDiffusionFlux(model, state, s, rho, yc, coeff1, Dij_row, poro, compIdx)
op = model.operators;
ncomp = model.getNumberOfComponents();
y = model.getProps(state, 'y');

% Millington-Quirk tortuosity
tau_mq = (s .* poro).^(7/3) .* poro.^(-2);

% Blanc's law for effective binary diffusion
D_diffij_inv = zeros(size(yc));
for cj = 1:ncomp
    if cj ~= compIdx
        if iscell(y)
            ycj = y{cj};
        else
            ycj = y(:, cj);
        end
        D_diffij_inv = D_diffij_inv + ycj ./ Dij_row(cj);
    end
end

D_diffij = coeff1 ./ max(D_diffij_inv, 1.e-16);
D_diff = op.faceAvg(s .* rho .* D_diffij .* tau_mq .* poro);

% Diffusion flux
flux = -D_diff .* op.Grad(yc);
end