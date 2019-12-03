classdef GenericBlackOilModel < ThreePhaseBlackOilModel & ExtendedReservoirModel
    properties

    end

    methods
        function model = GenericBlackOilModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});
            model.OutputStateFunctions = {'ComponentTotalMass'};
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@PhysicalModel(model, state0, state, dt, drivingForces, varargin{:});
        end

        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            [eqs, flux, names, types] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Treat source or bc terms
            if ~isempty(drivingForces.bc) || ~isempty(drivingForces.src)
                [pressures, sat, mob, rho, rs, rv] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'Rs', 'Rv');
                dissolved = model.getDissolutionMatrix(rs, rv);
                eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                                 pressures, sat, mob, rho, ...
                                                                 dissolved, {}, ...
                                                                 drivingForces);
            end
            % Assemble equations and add in sources
            for i = 1:numel(eqs)
                if ~isempty(src.cells)
                    eqs{i}(src.cells) = eqs{i}(src.cells) - src.value{i};
                end
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end
            % Get facility equations
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            eqs = [eqs, weqs];
            names = [names, wnames];
            types = [types, wtypes];
        end

        function names = getComponentNames(model)
            names = cellfun(@(x) x.name, model.Components, 'UniformOutput', false);
        end

        function [state, report] = updateState(model, state, problem, dx, forces)
            [state, report] = updateState@ThreePhaseBlackOilModel(model, state, problem, dx, forces);
            if ~isempty(model.FacilityModel)
                state = model.FacilityModel.applyWellLimits(state);
            end
        end

        function model = validateModel(model, varargin)
            % Validate model.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.validateModel`
            if isempty(model.FacilityModel) || ~isa(model.FacilityModel, 'ExtendedFacilityModel')
                model.FacilityModel = ExtendedFacilityModel(model);
            end
            if isempty(model.Components)
                nph = model.getNumberOfPhases();
                model.Components = cell(1, nph);
                names = model.getPhaseNames();
                disgas = model.disgas;
                vapoil = model.vapoil;
                for ph = 1:nph
                    switch names(ph)
                        case 'W'
                            c = ImmiscibleComponent('water', ph);
                        case 'O'
                            if disgas || vapoil
                                c = OilComponent('oil', ph, disgas, vapoil);
                            else
                                c = ImmiscibleComponent('oil', ph);
                            end
                        case 'G'
                            if disgas || vapoil
                                c = GasComponent('gas', ph, disgas, vapoil);
                            else
                                c = ImmiscibleComponent('gas', ph);
                            end
                        otherwise
                            error('Unknown phase');
                    end
                    model.Components{ph} = c;
                end
            end
            model = validateModel@ThreePhaseBlackOilModel(model, varargin{:});
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
            if model.outputFluxes
                state_flow = model.FluxDiscretization.buildFlowState(model, state, state0, dt);
                f = model.getProp(state_flow, 'PhaseFlux');
                nph = numel(f);
                state.flux = zeros(model.G.faces.num, nph);
                state.flux(model.operators.internalConn, :) = [f{:}];
            end
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
