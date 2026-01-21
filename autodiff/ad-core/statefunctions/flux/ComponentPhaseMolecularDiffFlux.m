classdef ComponentPhaseMolecularDiffFlux < StateFunction
    % ComponentPhaseMolecularDiffFlux - Molecular diffusion flux of each component in each phase
    %
    % Uses Quirk-Millington model for tortuosity with Blanc's law for gas
    % phase binary diffusion and literature values for liquid phase.
    %
    % Author: Stéphanie Delage Santacreu
    % Organization: Université de Pau et des Pays de l'Adour, E2S UPPA,
    %                CNRS, LFCR, UMR5150, Pau, France

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

            % Main calculation
            J = calculateMolecularDiffusionFluxes(prop, model, state);
        end
    end
end

%% Helper functions

function gp = setupDependencies(gp, model)
% Set up all required dependencies
gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
gp = gp.dependsOn('PoreVolume', 'PVTPropertyFunctions');
gp = gp.dependsOn('s', 'state');
gp = gp.dependsOn('x', 'state');
gp = gp.dependsOn('y', 'state');
gp = gp.dependsOn('pressure', 'state');
gp = gp.dependsOn('T', 'state');
end

function J = calculateMolecularDiffusionFluxes(prop, model, state)
% Calculate molecular diffusion fluxes for all components and phases

ncomp = model.getNumberOfComponents;
nph = model.getNumberOfPhases;
nm = model.getPhaseNames();

% Initialize output
J = cell(ncomp, nph);
[J{:}] = deal(0);

% Get diffusion coefficients
[mol_diff, param_LJ, Dij] = getDiffusionParameters(model);

% Get properties
rho = prop.getEvaluatedExternals(model, state, 'Density');
poro = model.getProp(state, 'Porosity');
[p, T] = model.getProps(state, 'pressure', 'temperature');
%  poro = getPorosity(model, p, state);

% Coefficients for gas diffusion
coeff1 = T.^(1.5) ./ (9.869e-6 .* p); % Convert pressure to atm

% Get phase indices
L_ix = model.getLiquidIndex();
V_ix = model.getVaporIndex();

% Calculate fluxes for each component and phase
for c = 1:ncomp
    [xc, yc] = getMoleFractions(state, c);

    for ph = 1:nph
        s = model.getProp(state, ['s', nm(ph)]);

        if ph == L_ix
            J{c, ph} = calculateLiquidDiffusionFlux(prop, model, s, rho{ph}, ...
                xc, mol_diff(c, ph), poro);

        elseif ph == V_ix
            J{c, ph} = calculateGasDiffusionFlux(prop, model, state, s, rho{ph}, ...
                yc, coeff1, Dij(c, :), poro, c);
        end
    end
end
end

function [mol_diff, param_LJ, Dij] = getDiffusionParameters(model)
% Load diffusion coefficients and calculate binary diffusion matrix

% Get component names
namecp = model.compFluid.names();

% Define component indices
componentNames = {'H2', 'C1', 'CO2', 'H2O', 'N2', 'C2', 'C3', 'H2S', 'NC4'};
indices = createComponentIndices(namecp, componentNames);

% Load diffusion coefficients (40°C)
[mol_diff, param_LJ] = loadDiffusionCoefficients(indices, componentNames);

% Calculate binary diffusion coefficients
Dij = calculateBinaryDiffusionMatrix(model, param_LJ);
end

function indices = createComponentIndices(namecp, componentNames)
% Create indices for component lookup
indices = struct();
for i = 1:numel(componentNames)
    name = componentNames{i};
    idx = find(strcmp(namecp, name));
    if ~isempty(idx)
        indices.(name) = idx;
    end
end
end

function [mol_diff, param_LJ] = loadDiffusionCoefficients(indices, componentNames)
% Load diffusion coefficients from databases

ncomp = max(structfun(@(x) x, indices));
mol_diff = zeros(ncomp, 2);
param_LJ = zeros(ncomp, 2);

% Diffusion coefficients at 40°C [liquid, gas] (m²/s)
% Note: Salinity reduces coefficients by 10-30%
coeffs = struct(...
    'H2',  [6.44e-9, 6.1e-5], ...
    'C1',  [2.15e-9, 1.6e-5], ...
    'H2O', [3.29e-9, 1.5e-5], ...
    'CO2', [2.72e-9, 1.4e-5], ...
    'N2',  [2.86e-9, 1.8e-5], ...
    'C2',  [1.72e-9, 2.5e-5], ...
    'C3',  [1.43e-9, 2.2e-5], ...
    'H2S', [2.15e-9, 2.2e-5], ...
    'NC4', [1.15e-9, 1.9e-5]);

% Lennard-Jones parameters [diameter (Å), potential (K)]
coeffs_LJ = struct(...
    'H2',  [2.92,   59.7], ...
    'C1',  [3.758, 148.6], ...
    'H2O', [2.641, 809.1], ...
    'CO2', [3.996, 195.2], ...
    'N2',  [3.798,  71.4], ...
    'C2',  [4.443, 215.7], ...
    'C3',  [5.118, 237.1], ...
    'H2S', [3.60,  301.0], ...
    'NC4', [5.206, 289.5]);

% Assign coefficients to components
for i = 1:numel(componentNames)
    comp = componentNames{i};
    if isfield(indices, comp) && isfield(coeffs, comp)
        idx = indices.(comp);
        mol_diff(idx, :) = coeffs.(comp);
        param_LJ(idx, :) = coeffs_LJ.(comp);
    end
end
end

function Dij = calculateBinaryDiffusionMatrix(model, param_LJ)
% Calculate binary diffusion coefficients matrix

ncomp = size(param_LJ, 1);
Molmass = 1.e3 .* model.compFluid.molarMass;
SigLJ = param_LJ(:, 1);

Dij = zeros(ncomp);
for c = 1:ncomp
    for cj = 1:ncomp
        sqrtMij = sqrt(2 * Molmass(c) * Molmass(cj) / (Molmass(c) + Molmass(cj)));
        Sigij2 = 0.25 * (SigLJ(c) + SigLJ(cj))^2;
        Dij(c, cj) = 1.e-4 * 0.001858 / (sqrtMij * Sigij2); % m²/s
    end
end
end

function poro = getPorosity(model, p, state)
% Get porosity with optional bioclogging
if model.dynamicFlowPv()
    if model.bacteriamodel
        poro = model.rock.poro(p, state.nbact);
    else
        poro = model.rock.poro(p);
    end
else
    poro = model.rock.poro;
end
end

function [xc, yc] = getMoleFractions(state, c)
% Extract mole fractions for component c
if iscell(state.x)
    xc = state.x{c};
else
    xc = state.x(:, c);
end

if iscell(state.y)
    yc = state.y{c};
else
    yc = state.y(:, c);
end
end

function flux = calculateLiquidDiffusionFlux(prop, model, s, rho, xc, D_mol, poro)
% Calculate liquid phase diffusion flux
avg = model.operators.faceAvg;

% Millington-Quirk tortuosity model
tau_mq = (s .* poro).^(7/3) .* poro.^(-2);

% Effective diffusivity
D_eff = avg(s .* rho .* D_mol .* tau_mq .* poro);

% Diffusion flux
flux = -D_eff .* model.operators.Grad(xc);
end

function flux = calculateGasDiffusionFlux(prop, model, state, s, rho, yc, coeff1, Dij_row, poro, compIdx)
% Calculate gas phase diffusion flux using Blanc's law
avg = model.operators.faceAvg;
ncomp = model.getNumberOfComponents();

% Millington-Quirk tortuosity model
tau_mq = (s .* poro).^(7/3) .* poro.^(-2);

% Calculate effective binary diffusion coefficient
D_diffij_inv = zeros(size(yc));
for cj = 1:ncomp
    if cj ~= compIdx
        if iscell(state.y)
            ycj = state.y{cj};
        else
            ycj = state.y(:, cj);
        end
        D_diffij_inv = D_diffij_inv + ycj ./ Dij_row(cj);
    end
end

D_diffij = coeff1 ./ max(D_diffij_inv, 1.e-16);
D_diff = avg(s .* rho .* D_diffij .* tau_mq .* poro);

% Diffusion flux
flux = -D_diff .* model.operators.Grad(yc);
end

%% Coefficient tables (for reference)
%{
% Diffusion coefficients at 25°C (commented out for reference)
% coeffs_25C = struct(...
%     'H2',  [4.5e-9, 6.1e-5], ...
%     'C1',  [1.5e-9, 1.6e-5], ...
%     'H2O', [2.3e-9, 1.5e-5], ...
%     'CO2', [1.9e-9, 1.4e-5], ...
%     'N2',  [2.0e-9, 1.8e-5], ...
%     'C2',  [1.2e-9, 2.5e-5], ...
%     'C3',  [1.0e-9, 2.2e-5], ...
%     'H2S', [1.5e-9, 2.2e-5], ...
%     'NC4', [0.8e-9, 1.9e-5]);
%}

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}