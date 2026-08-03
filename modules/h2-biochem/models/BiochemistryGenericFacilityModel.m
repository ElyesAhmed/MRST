classdef BiochemistryGenericFacilityModel < GenericFacilityModel
    % BiochemistryGenericFacilityModel
    % Generic facility model for biochemistry simulations, including
    % bacterial growth and decay.

    properties
        bacterialFormulation = 'bacterialmodel';  % Formulation for bacterial transport
    end

    methods
        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, useDefaults)
            % Set up state function groupings using parent, and add
            % biochemistry-specific functions

            if nargin < 2
                useDefaults = isempty(model.FacilityFlowDiscretization);
            end

            % Base facility groupings
            model = setupStateFunctionGroupings@GenericFacilityModel(model, useDefaults);

            % Add biochemistry-specific state functions only when bacteriamodel is active
            rm = model.ReservoirModel;
            if ~isempty(rm) && isprop(rm, 'bacteriamodel') && rm.bacteriamodel
                ffd = model.FacilityFlowDiscretization;
                % Note: BacterialMass is already registered in ReservoirModel's PVTPropertyFunctions
                ffd = ffd.setStateFunction('PsiGrowthRate', GrowthBactRateSRC(model));
                ffd = ffd.setStateFunction('PsiDecayRate', DecayBactRateSRC(model));
                ffd = ffd.setStateFunction('BactConvRate', BactConvertionRate(model));
                if isprop(rm, 'sulfateReduction') && rm.sulfateReduction
                    ffd = ffd.setStateFunction('SRBTracerConvRate', SRBTracerConvRate(model));
                end
                model.FacilityFlowDiscretization = ffd;
            end
        end

        %-----------------------------------------------------------------%
        function state = initStateAD(model, state, vars, names, origin)
            % Initialize AD state from double state
            state = initStateAD@GenericFacilityModel(model, state, vars, names, origin);
        end

        %-----------------------------------------------------------------%
        function names = getBasicPrimaryVariableNames(model)
            % Get names of primary variables
            names = getBasicPrimaryVariableNames@GenericFacilityModel(model);

            % If bacterial variables are disabled, return immediately
            if strcmpi(model.bacterialFormulation, 'none') || ...
                    strcmpi(model.primaryVariableSet, 'none')
                return
            end
        end

        %-----------------------------------------------------------------%
        function [variables, names, map] = getBasicPrimaryVariables(model, wellSol)
            % Return facility primary variables
            [variables, names, map] = getBasicPrimaryVariables@GenericFacilityModel(model, wellSol);
        end

        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            % Get field name and index for a given variable
            [fn, index] = getVariableField@GenericFacilityModel(model, name, varargin{:});
        end

        %-----------------------------------------------------------------%
        function src = getBacteriaSources(model, fd, state, state0, dt)
            % Growth/decay source: src = (g - d) * BacterialMass
            % where g, d are kinetic rates (1/s) applied to mass directly
            rm = model.ReservoirModel;
            bcrm=rm.biochemFluid;
            if isempty(rm) || ~isprop(rm, 'bacteriamodel') || ~rm.bacteriamodel
                src = 0;
                return;
            end

            reg = 1.0e-10;
            flowState = fd.buildFlowState(model, state, state0, dt);
            psigrowth = model.getProps(flowState, 'PsiGrowthRate');  % Psigrowthmax * axH2 * axsub [1/s]
            psidecay  = model.getProps(flowState, 'PsiDecayRate');   % bbact * nbact [1/s]
            bmass     = rm.PVTPropertyFunctions.get(rm, state, 'BacterialMass');  % pv * S_l * rho_l * nbact [kg]

            % Direct (g-d)*mass formulation using BacterialMass
            nbioreact=bcrm.nbioreact;
            src_growthdecay = cell(1,nbioreact);
            [src_growthdecay{:}] = deal(0);
            for i=1:nbioreact
                src_growthdecay{i} = (psigrowth{i} - psidecay{i}).* bmass{i} - reg .* bmass{i};
            end

            % ===== NEW: Well bacteria source (advective transport) =====
            map   = model.getProp(state, 'FacilityWellMapping');
            src_well = cell(1,nbioreact);
            [src_well{:}] = deal(0);
            if ~isempty(map.cells)
                q_ph  = model.getProp(state, 'PhaseFlux');
                rho = rm.PVTPropertyFunctions.get(rm, state, 'Density');
                nbact = rm.getProp(state, 'nbact');
                L_ix  = rm.getLiquidIndex();

                for i=1:nbioreact
                    % Liquid phase flux per perforation (positive = injection)
                    q_l   = q_ph{L_ix};
                    % Liquid density in perforated cells
                    rho_l = rho{L_ix};
                    % Bacteria mass flux: ρ_l * q_l * ω
                    rho_perf = rho_l(map.cells);
                    if iscell(nbact)
                        nbacti=nbact{i};
                    else
                        nbacti=nbact(:,i);
                    end
                    q_bact = rho_perf .* q_l .* nbacti(map.cells);
                    % Injectors: no bacteria injected → set to 0
                    q_bact(q_l > 0) = 0;

                    % Sum perforation contributions to cells (producers give negative)
                    % Sparse summation to cells (AD‑compatible)
                    nc = rm.G.cells.num;
                    S = sparse(map.cells, (1:numel(map.cells))', 1, nc, numel(map.cells));
                    src_well{i} = S * q_bact;
                end

            end
            % ==============================================================
            src = cell(1,nbioreact);
            for i=1:nbioreact
                src{i} = src_growthdecay{i} + src_well{i};
            end
        end
        %-----------------------------------------------------------------%
        function src = getSRBTracerSources(model, fd, state, state0, dt)
            % Reaction + well-advection source for the SO4/HS aqueous
            % tracers. Mirrors getBacteriaSources, but for the two
            % non-volatile species that never enter the EOS/flash.
            rm = model.ReservoirModel;
            if isempty(rm) || ~isprop(rm, 'sulfateReduction') || ~rm.sulfateReduction
                src = {0, 0};
                return;
            end

            flowState = fd.buildFlowState(model, state, state0, dt);
            src_reaction = model.getProps(flowState, 'SRBTracerConvRate');

            % Well tracer source (advective transport, no SO4/HS injected)
            map = model.getProp(state, 'FacilityWellMapping');
            src_well = {0, 0};
            if ~isempty(map.cells)
                q_ph = model.getProp(state, 'PhaseFlux');
                L_ix = rm.getLiquidIndex();
                q_l  = q_ph{L_ix};

                so4 = rm.getProp(state, 'so4');
                hs  = rm.getProp(state, 'hs');
                tracers = {so4, hs};

                nc = rm.G.cells.num;
                nf = numel(map.cells);
                S  = sparse(map.cells, (1:nf)', 1, nc, nf);

                for i = 1:2
                    q_tracer = q_l .* tracers{i}(map.cells);
                    % Injectors: no SO4/HS injected -> set to 0
                    q_tracer(q_l > 0) = 0;
                    src_well{i} = S * q_tracer;
                end
            end

            src = cell(1, 2);
            for i = 1:2
                src{i} = src_reaction{i} + src_well{i};
            end
        end

        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            % Return facility equations including parent contributions
            [eqs, names, types, state] = ...
                getModelEquations@GenericFacilityModel(model, state0, state, dt, drivingForces);
        end

        %-----------------------------------------------------------------%
        function [values, tolerances, names, evaluated] = getFacilityConvergenceValues(model, problem, varargin)
            % Get convergence values for facility
            [values, tolerances, names, evaluated] = ...
                getFacilityConvergenceValues@GenericFacilityModel(model, problem, varargin{:});
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