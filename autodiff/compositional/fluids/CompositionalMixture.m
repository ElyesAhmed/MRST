classdef CompositionalMixture
    properties
        Tcrit % Critical temperature in Kelvin
        Pcrit % Critical pressure in Pascal
        Vcrit % Critical volume in m^3 / mol
        acentricFactors % Acentric factors (dimensionless)
        molarMass % Component mass (kg / mol)
        names % Names of each component. Each name must be unique.
        name % Name of the fluid mixture itself
    end
    
    properties ( Access = protected)
        bic % Binary interaction coefficients
    end
    
    methods
        function fluid = CompositionalMixture(names, Tcrit, Pcrit, Vcrit, acentricFactors, molarMass, varargin)
            if nargin == 0
                disp('Empty fluid - no validation');
                return
            end
            
            for i = 1:numel(names)
                cts = sum(strcmp(names{i}, names));
                if cts > 1
                    warning(['Component ', num2str(i), ': ', names{i}, ' occurs multiple times.']);
                end
            end
            fluid.names = names;
            fluid.Tcrit = Tcrit;
            fluid.Pcrit = Pcrit;
            fluid.Vcrit = Vcrit;
            fluid.acentricFactors = acentricFactors;
            fluid.molarMass = molarMass;
            ncomp = numel(fluid.names);
            
            fluid.bic = zeros(ncomp, ncomp);
        end
        
        function n = getNumberOfComponents(fluid)
            n = numel(fluid.names);
        end
        
        function bic = getBinaryInteraction(fluid)
            bic = fluid.bic;
        end

%         function bic = getBinaryInteractionGH2H2O(fluid, T)
%             % Calculates gas-phase binary interaction coefficients (BIC)
%             % between H2 and H2O.
% 
%             % Extract component indices
%             namecp = fluid.names;
%             indH2 = find(strcmp(namecp, 'H2'));
%             indH2O = find(strcmp(namecp, 'H2O'));
% 
%             % Compute reduced temperature for H2
%             TrH2 = T ./ fluid.Tcrit(indH2);
% 
%             % Gas-phase BIC coefficients for H2-H2O
%             [D0, D1] = deal(0.01993, 0.042834);
%             bic{indH2, indH2O} = D0 + D1 .* TrH2;
%             bic{indH2O, indH2} = bic{indH2, indH2O}; % Symmetry
% 
%             % Self-interactions (diagonal terms)
%             bic{indH2, indH2} = 0 .* TrH2 + fluid.bic(indH2, indH2);
%             bic{indH2O, indH2O} = 0 .* TrH2 + fluid.bic(indH2O, indH2O);
%         end
% 
%         function bic = getBinaryInteractionLH2H2O(fluid, T, msalt)
%             % Calculates liquid-phase binary interaction coefficients (BIC)
%             % between H2 and H2O with salinity effects.
% 
%             % Extract component indices
%             namecp = fluid.names;
%             indH2 = find(strcmp(namecp, 'H2'));
%             indH2O = find(strcmp(namecp, 'H2O'));
% 
%             % Compute reduced temperature for H2
%             TrH2 = T ./ fluid.Tcrit(indH2);
% 
%             % Liquid-phase BIC coefficients for H2-H2O with salinity influence
%             [D0, D1, D2, D3] = deal(-2.11917, 0.14888, -13.01835, -0.43946);
%             [a0, a1] = deal(-2.26322e-2, -4.4736e-3);
% 
%             bic{indH2, indH2O} = (D0 .* (1 + a0 .* msalt) + ...
%                 D1 .* TrH2 .* (1 + a1 .* msalt) + ...
%                 D2 .* exp(D3 .* TrH2));
%             bic{indH2O, indH2} = bic{indH2, indH2O}; % Symmetry
% 
%             % Self-interactions (diagonal terms)
%             bic{indH2, indH2} = 0 .* TrH2 + fluid.bic(indH2, indH2);
%             bic{indH2O, indH2O} = 0 .* TrH2 + fluid.bic(indH2O, indH2O);
%         end


function bic = getBinaryInteractionGasWater(fluid, T)
    % Generalized gas-phase BIC calculation for components with H2O.
    %
    % Inputs:
    %   - fluid: structure containing component names, critical temperatures, and BICs
    %   - T: temperature (scalar or vector)
    %
    % Outputs:
    %   - bic: binary interaction coefficients matrix

    % Extract component indices dynamically
    namecp = fluid.names;
    indices = struct('H2', find(strcmp(namecp, 'H2')), ...
                     'CH4', find(strcmp(namecp, 'CH4')), ...
                     'CO2', find(strcmp(namecp, 'CO2')), ...
                     'N2', find(strcmp(namecp, 'N2')), ...
                     'H2O', find(strcmp(namecp, 'H2O')));

    % Initialize BIC matrix
    nComponents = numel(namecp);
    bic = cell(nComponents);

    % Define gas-phase BIC coefficients
    coeffs = struct('H2', [0.01993, 0.042834], ...
        'CH4', [0.025, 0.035], ...
        'CO2', [0.018, 0.028], ...
        'N2', [0.022, 0.032]);

    % Assign BIC values for each component interacting with H2O
    fields = fieldnames(indices);
    for i = 1:numel(fields)
        comp = fields{i};
        indComp = indices.(comp);
        indH2O = indices.H2O;

        if ~isempty(indComp) && ~isempty(indH2O) && isfield(coeffs, comp)
            Tr = T ./ fluid.Tcrit(indComp);
            bic{indComp, indH2O} = coeffs.(comp)(1) + coeffs.(comp)(2) .* Tr;
            bic{indH2O, indComp} = bic{indComp, indH2O}; % Symmetry
        end
    end
    % Handle self-interactions for H2, CO2, CH4, N2, and H2O
    for i=1:nComponents
        bic{i,i} = fluid.bic(i,i)+0.*Tr;
    end
end
function bic = getBinaryInteractionLiquidWater(fluid, T, msalt)
    % Generalized liquid-phase BIC calculation for components with H2O.
    %
    % Inputs:
    %   - fluid: structure containing component names, critical temperatures, and BICs
    %   - T: temperature (scalar or vector)
    %   - msalt: salinity (scalar or vector)
    %
    % Outputs:
    %   - bic: binary interaction coefficients matrix

    % Extract component indices dynamically
    namecp = fluid.names;
    indices = struct('H2', find(strcmp(namecp, 'H2')), ...
                     'CH4', find(strcmp(namecp, 'CH4')), ...
                     'CO2', find(strcmp(namecp, 'CO2')), ...
                     'N2', find(strcmp(namecp, 'N2')), ...
                     'H2O', find(strcmp(namecp, 'H2O')));

    % Initialize BIC matrix
    nComponents = numel(namecp);
    bic = cell(nComponents);

    % Define liquid-phase BIC coefficients with salinity influence
    coeffs = struct('H2', [-2.11917, 0.14888, -13.01835, -0.43946, -2.26322e-2, -4.4736e-3], ...
                    'CH4', [-1.5, 0.12, -10, -0.35, -0.02, -0.004], ...
                    'CO2', [-1.8, 0.13, -11.5, -0.4, -0.025, -0.005], ...
                    'N2', [-2.0, 0.15, -12, -0.45, -0.03, -0.006]);

    % Assign BIC values for each component interacting with H2O
    fields = fieldnames(indices);
    for i = 1:numel(fields)
        comp = fields{i};
        indComp = indices.(comp);
        indH2O = indices.H2O;

        if ~isempty(indComp) && ~isempty(indH2O) && isfield(coeffs, comp)
            Tr = T ./ fluid.Tcrit(indComp);
            bic{indComp, indH2O} = coeffs.(comp)(1) .* (1 + coeffs.(comp)(5) .* msalt) + ...
                                   coeffs.(comp)(2) .* Tr .* (1 + coeffs.(comp)(6) .* msalt) + ...
                                   coeffs.(comp)(3) .* exp(coeffs.(comp)(4) .* Tr);
            bic{indH2O, indComp} = bic{indComp, indH2O}; % Symmetry
        end
    end
    % Handle self-interactions for H2, CO2, CH4, N2, and H2O
    for i=1:nComponents
        bic{i,i} = fluid.bic(i,i)+0.*Tr;
    end
end


        function fluid = setBinaryInteraction(fluid, input)
            % Set BIC via a matrix. Must be symmetric and ncomp by ncomp
            ncomp = fluid.getNumberOfComponents();
            
            if isvector(input)
                n_el = ncomp*(ncomp-1)/2;
                assert(numel(input) == n_el);
                
                tmp = zeros(ncomp, ncomp);
                offset = 0;
                for i = 2:ncomp
                    cts = i-1;
                    tmp(i, 1:i-1) = input(offset + (1:cts));
                    offset = offset + cts;
                end
                input = tmp + tmp';
            end
            assert(size(input, 1) == ncomp)
            assert(size(input, 1) == size(input, 2));
            assert(all(all(input == input')));
            fluid.bic = input;
        end
        
        function disp(mixture)
            % builtin('disp', mixture);
            cn = class(mixture);
            isDesktop = usejava('desktop');
            if isDesktop
                fprintf('  <a href="matlab:helpPopup %s">%s</a>:\n', cn, cn);
            else
                fprintf('  %s:\n', cn);
            end
            fprintf('\n');
            nc = mixture.getNumberOfComponents();
            cnames = mixture.names;
            fprintf('  %d component mixture', nc);
            if ~isempty(mixture.name)
                fprintf(' (%s)', mixture.name)
            end
            fprintf(':\n');

            ml = max(max(cellfun(@numel, cnames)), 4);
            sep = repmat('-', 55 + ml, 1);
            fprintf('  %*s | p_c [Pa] | T_c [K] | V_c [m^3] |  acf  | mw [kg/mol] \n', ml, 'Name');
            fprintf('  %s\n', sep)
            for i = 1:nc
                pc = mixture.Pcrit(i);
                tc = mixture.Tcrit(i);
                vc = mixture.Vcrit(i);
                acf = mixture.acentricFactors(i);
                mw = mixture.molarMass(i);
                fprintf('  %*s | %1.2e | %3.1f K | %2.3e | %1.3f | %1.7f \n', ...
                        ml, cnames{i}, pc, tc, vc, acf, mw);
            end
            fprintf('  %s\n', sep);
            bi = mixture.getBinaryInteraction();
            fprintf('  ');
            if any(bi(:))
                fprintf('Binary interaction coefficients:\n');
                disp(bi);
            else
                fprintf('No non-zero binary interaction coefficients.\n')
            end
            fprintf('\n');
            known = {'Pcrit', 'Tcrit', 'Vcrit', 'acentricFactors', 'molarMass', 'names', 'name'};
            props = propertynames(mixture);
            extra = setdiff(props, known);
            if numel(extra)
                fprintf('\n  Additional properties:\n');
                for i = 1:numel(extra)
                    e = extra{i};
                    fprintf('  %s (%s)\n', e, class(mixture.(e)));
                end
                fprintf('\n');
            end
        end
    end
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
