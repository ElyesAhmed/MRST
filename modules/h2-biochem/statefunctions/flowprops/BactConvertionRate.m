classdef BactConvertionRate < StateFunction
    % Bacterial conversion rate model for compositional simulations
    %
    % SYNOPSIS:
    %   gf = BactConvertionRate(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   This class computes the bacterial conversion rate for components in a
    %   compositional simulation with microbial activity. The model accounts
    %   for bacterial growth limited by H2 and CO2 availability and converts
    %   these substrates into biomass and products according to stoichiometric
    %   coefficients.
    %
    % SEE ALSO:
    %   CompositionalModel, EquationsCompositional

    properties
        % No additional properties needed - all parameters come from model
    end

    methods
        function gp = BactConvertionRate(model, varargin)
            % Constructor for bacterial conversion rate calculator
            gp@StateFunction(model, varargin{:});

            % Define dependencies
            gp = gp.dependsOn({'PsiGrowthRate', 'PsiDecayRate'}, 'FlowDiscretization');
            gp = gp.dependsOn({'PoreVolume'}, 'PVTPropertyFunctions');

            % Set label for output
            gp.label = 'Q_biot'; %
        end

        function qbiot = evaluateOnDomain(prop, model, state)
            % Compute bacterial conversion rate for each component
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Compositional model instance
            %   state - State struct with fields
            %
            % RETURNS:
            %   qbiot - Cell array of conversion rates per component [kg/s]

            % Initialize with zeros
            ncomp = model.ReservoirModel.EOSModel.getNumberOfComponents();
            qbiot = cell(ncomp, 1);
            [qbiot{:}] = deal(0);

            % Check if bacterial modeling is active
            rm = model.ReservoirModel;
            if ~(rm.bacteriamodel && rm.liquidPhase)
                return;
            end
            bcrm=rm.biochemFluid;

            % Get component indices
            compNames = rm.EOSModel.getComponentNames();
            idxH2  = find(strcmpi(compNames, bcrm.rH2));
            idxsub = find(strcmpi(compNames, bcrm.rsub));
            

            % Validate required components
            if strcmp(bcrm.metabolicReaction, 'MethanogenicArchae')
                idxCO2 = find(strcmpi(compNames, 'CO2'), 1);
                if isempty(idxH2) || isempty(idxCO2)
                    warning('Bacterial model requires H2 and CO2 components');
                    return;
                end
            end
            if strcmp(bcrm.metabolicReaction, 'SulfateReducingBacteria')
                idxSO4 = find(strcmpi(compNames, 'SO4'), 1);
                if isempty(idxH2) || isempty(idxSO4)
                    warning('Bacterial model requires H2 and SO4(2-) components');
                    return;
                end
            end
            try
                % Get required state variables
                pv = rm.PVTPropertyFunctions.get(model.ReservoirModel, state, 'PoreVolume');
                s  =  rm.getProp(state, 's');
                nbact = rm.getProp(state, 'nbact');
                L_ix = rm.getLiquidIndex();
                x = rm.getProp(state, 'x');

                % Extract liquid phase properties
                if iscell(x)
                    xH2 = x{idxH2};
                    xsub = x{idxsub};
                    sL = s{L_ix};
                else
                    xH2 = x(:, idxH2);
                    xsub = x(:, idxsub);
                    sL = s(:, L_ix);
                end

                % Ensure non-zero liquid saturation
                sL = max(value(sL), 1.0e-8);

                % Get model parameters
                alphaH2 = bcrm.alphaH2;
                alphasub = bcrm.alphasub;
                Psigrowthmax = bcrm.Psigrowthmax;
                Y_H2 = bcrm.Y_H2;
                gammak = rm.gammak;
                nbactMax = bcrm.nbactMax;
                mc = rm.EOSModel.CompositionalMixture.molarMass;

                % Calculate growth rate using Monod kinetics
                axH2 = xH2 ./ (alphaH2 + xH2);
                axsub = xsub ./ (alphasub + xsub);

                Psigrowth = pv .* Psigrowthmax .* axH2 .* axsub .* nbact .* sL;

                % Calculate conversion rates for all components
                qbiot_temp = Psigrowth ./ Y_H2 ./ abs(gammak(idxH2));
                for c = 1:ncomp
                    qbiot{c} = gammak(c) .* qbiot_temp .* mc(c) .* nbactMax;
                end

            catch ME
                warning('Bacterial conversion rate calculation failed: %s', ME.message);
                [qbiot{:}] = deal(0);
            end
        end
    end
end
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