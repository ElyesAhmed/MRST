classdef ComponentPhaseMolecularDiffFlux < StateFunction
    properties
        D                        % Face diffusivity
        diffusionCoefficients    % Component diffusion coefficients [liquid, gas]
        tortuosityFactor = 7/3   % Millington-Quirk exponent
        useFaceAverage = true    % Use face average (false for upwinding)
    end

    methods
        function cf = ComponentPhaseMolecularDiffFlux(model)
            cf@StateFunction(model);

            % Set up dependencies
            cf = cf.dependsOn({'Density', 's', 'pressure'}, 'state');
            if model.isCompositional
                cf = cf.dependsOn({'x', 'y'}, 'state');
            else
                cf = cf.dependsOn({'rs', 'rv', 'SurfaceDensity'}, 'state');
            end

            % Pre-compute face diffusivity
            cf.D = getFaceDiffusivity(model.G, model.rock);

            % Get component names and setup diffusion coefficients
            ncomp = model.getNumberOfComponents();
            componentNames = cell(ncomp, 1);
            for c = 1:ncomp
                componentNames{c} = model.Components{c}.name;
            end

            % Load diffusion coefficients
            cf.diffusionCoefficients = loadDiffusionCoefficients(componentNames);
        end

        function v = evaluateOnDomain(cf, model, state)
            % Get the convective flux from parent class
            v = evaluateOnDomain@StateFunction(cf, model, state);

            op = model.operators;
            T = cf.D(op.internalConn);

            % Get phases
            nph = model.getNumberOfPhases();
            L_ix = model.getLiquidIndex();
            V_ix = model.getVaporIndex();

            % Get required properties
            rho = model.getProps(state, 'Density');
            s = state.s;
            p = state.pressure;

            % Get porosity
            if isfield(model.rock, 'poro')
                poro = model.rock.poro;
            else
                poro = model.rock.poro * ones(model.G.cells.num, 1);
            end

            % Get mole fractions (x for liquid, y for gas)
            if model.isCompositional
                x = state.x;
                y = state.y;
            else
                % For black-oil, calculate from rs/rv
                [x, y] = calculateBlackOilMoleFractions(model, state);
            end

            % Add diffusion fluxes for each component
            ncomp = model.getNumberOfComponents();
            for c = 1:ncomp
                % Get appropriate mole fraction for each phase
                if iscell(x)
                    xc = x{c};
                    yc = y{c};
                else
                    xc = x(:, c);
                    yc = y(:, c);
                end

                for ph = 1:nph
                    if ph == L_ix
                        % Liquid phase diffusion
                        flux = calculateLiquidDiffusionFlux(cf, model, state, ...
                            c, ph, xc, rho{ph}, s{ph}, poro, T);

                    elseif ph == V_ix
                        % Gas phase diffusion
                        flux = calculateGasDiffusionFlux(cf, model, state, ...
                            c, ph, yc, y, rho{ph}, s{ph}, poro, p, state.temperature, T);
                    end

                    % Add diffusion flux to convective flux
                    v{c, ph} = v{c, ph} + flux;
                end
            end
        end
    end
end

%% Helper functions (simplified)

function [x, y] = calculateBlackOilMoleFractions(model, state)
% Simplified version for black-oil
ncomp = model.getNumberOfComponents();
x = cell(ncomp, 1);
y = cell(ncomp, 1);

% Get rs and rv
rs = getValue(state, 'rs', 0);
rv = getValue(state, 'rv', 0);

% Simple approximation
for c = 1:ncomp
    % In oil phase
    x{c} = 1 ./ (1 + rs);

    % In gas phase
    y{c} = rv ./ (1 + rv);
end

% Adjust for oil/gas components
if ncomp >= 2
    x{2} = rs ./ (1 + rs);  % Gas in oil
    y{1} = rv ./ (1 + rv);  % Oil in gas
end
end

function flux = calculateLiquidDiffusionFlux(cf, model, state, compIdx, phIdx, xc, rho, s, poro, T)
op = model.operators;

% Tortuosity factor
tau = (s .* poro).^(cf.tortuosityFactor) .* poro.^(-2);

% Get diffusion coefficient
D_mol = cf.diffusionCoefficients(compIdx, 1);

% Calculate gradient
grad_x = op.Grad(xc);

if cf.useFaceAverage
    faceDens = op.faceAvg(rho);
    faceS = op.faceAvg(s);
    faceTau = op.faceAvg(tau);
    Tc = D_mol .* T .* faceDens .* faceS .* faceTau;
else
    % Upwinding based on gradient direction
    flag = value(grad_x) < 0;
    faceDens = op.faceUpstr(flag, rho);
    faceS = op.faceUpstr(flag, s);
    faceTau = op.faceUpstr(flag, tau);
    Tc = D_mol .* T .* faceDens .* faceS .* faceTau;
end

flux = -Tc .* grad_x;
end

function val = getValue(state, fieldname, default)
if isfield(state, fieldname)
    val = state.(fieldname);
else
    val = default;
end
end

%% Coefficient tables

function coeffs = loadDiffusionCoefficients(componentNames)
% Diffusion coefficients at 40°C [liquid, gas] (m²/s)
database = struct(...
    'H2',  [6.44e-9, 6.1e-5], ...
    'C1',  [2.15e-9, 1.6e-5], ...
    'H2O', [3.29e-9, 1.5e-5], ...
    'CO2', [2.72e-9, 1.4e-5], ...
    'N2',  [2.86e-9, 1.8e-5], ...
    'C2',  [1.72e-9, 2.5e-5], ...
    'C3',  [1.43e-9, 2.2e-5], ...
    'H2S', [2.15e-9, 2.2e-5], ...
    'NC4', [1.15e-9, 1.9e-5], ...
    'oil', [0.8e-9, 1.2e-5], ...
    'gas', [1.5e-9, 1.8e-5], ...
    'water', [2.3e-9, 1.5e-5]);

ncomp = numel(componentNames);
coeffs = zeros(ncomp, 2);

for c = 1:ncomp
    name = componentNames{c};
    if isfield(database, name)
        coeffs(c, :) = database.(name);
    else
        % Default values for unknown components
        if contains(lower(name), 'oil')
            coeffs(c, :) = [0.8e-9, 1.2e-5];
        elseif contains(lower(name), 'gas')
            coeffs(c, :) = [1.5e-9, 1.8e-5];
        elseif contains(lower(name), 'water')
            coeffs(c, :) = [2.3e-9, 1.5e-5];
        else
            coeffs(c, :) = [1e-9, 1e-5];
        end
    end
end
end
%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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
