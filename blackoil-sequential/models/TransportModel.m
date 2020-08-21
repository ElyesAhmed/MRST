classdef TransportModel < WrapperModel
    properties
        formulation = 'totalSaturation'
        dsMaxTotal = inf;
    end
    
    methods
        function model = TransportModel(parent, varargin)
            parent.FlowDiscretization = []; % Remove flux discretization (if setup)
            if isprop(parent, 'FacilityModel') && ~isempty(parent.FacilityModel)
                parent.FacilityModel.primaryVariableSet = 'none';
            end
            model = model@WrapperModel(parent);
            model = merge_options(model, varargin{:});
            model.AutoDiffBackend = parent.AutoDiffBackend;
        end
        
        function [state, names, origin] = getStateAD(model, state, init)
            if nargin < 3
                init = true;
            end
            parent = model.parentModel;
            % Get the AD state for this model
            [basevars, basenames, baseorigin] = model.getPrimaryVariables(state);
            isParent = strcmp(baseorigin, class(parent));
            isPrimaryVariable = isParent;
            % Find saturations
            isS = false(size(basevars));
            nph = parent.getNumberOfPhases();
            phase_variable_index = zeros(nph, 1);
            for i = 1:numel(basevars)
                if ~isParent(i)
                    % Skip wells etc
                    continue
                end
                bn = basenames{i};
                if strcmp(bn, 'x')
                    % Hack for blackoil-models
                    ss = true;
                else
                    [f, ix] = model.getVariableField(basenames{i});
                    ss = strcmp(f, 's');
                end
                if ss
                    isS(i) = true;
                    phase_variable_index(ix) = i;
                end
            end
            % Figure out saturation logic
            isP = strcmp(basenames, 'pressure');
            vars = basevars;
            names = basenames;
            origin = baseorigin;
            useTotalSaturation = strcmpi(model.formulation, 'totalSaturation') ...
                                     && (sum(isS) < nph);
            if useTotalSaturation
                % Replace pressure with total saturation
                replacement = 'sT';
                sT = model.getProp(state, replacement);
                % Replacing
                vars{isP} = sT; % {'pressure', 'sw', 'x'} -> {'st', 'sw', 'x'}
                names{isP} = replacement;
                origin{isP} = class(model);
            else
                % Remove pressure and skip saturation closure
                keep = ~isP;
                vars = vars(keep);
                names = names(keep);
                origin = origin(keep);
                isPrimaryVariable = isPrimaryVariable(keep);
            end
            if init
                [vars{isPrimaryVariable}] = model.AutoDiffBackend.initVariablesAD(vars{isPrimaryVariable});
            end
            if useTotalSaturation
                basevars(~isP) = vars(~isP);
            else
                basevars(~isP) = vars;
            end
            state = model.initStateAD(state, basevars, basenames, baseorigin);
            if useTotalSaturation
                % Set total saturation as well
                sT = vars{isP};
                state = model.setProp(state, replacement, sT);
            end
            % Finally remove the non-primary values from the list
            names = names(isPrimaryVariable);
            origin = origin(isPrimaryVariable);
        end
        

        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            parent = model.parentModel;
            [eqs, names, types, state] = parent.getModelEquations(state0, state, dt, drivingForces);
            if strcmpi(model.formulation, 'missingPhase')
                % Skip the last pseudocomponent/phase! Only
                % mass-conservative for incompressible problems or with
                % outer loop enabled
                cnames = parent.getComponentNames();
                subs = true(size(names));
                subs(strcmpi(names, cnames{end})) = false;
                eqs = eqs(subs);
                names = names(subs);
                types = types(subs);
            end
        end

        function state = validateState(model, state)
            state = validateState@WrapperModel(model, state);
            fn = model.getVariableField('sT');
            if isfield(state, 's') && ~isfield(state, fn)
                state = model.setProp(state, 'sT', sum(state.s, 2));
            end
        end
        
        function model = setupStateFunctionGroupings(model, setDefaults)
            if nargin < 2
                setDefaults = true;
            end
            pmodel = model.parentModel;
            hasFacility = isprop(pmodel, 'FacilityModel') && ~isempty(pmodel.FacilityModel);
            isTotalSat = strcmpi(model.formulation, 'totalSaturation');
            if isempty(pmodel.FlowDiscretization)
                pmodel = pmodel.setupStateFunctionGroupings();
            elseif ~setDefaults
                % State functions have already been set up
                warning('State function models for parent model is already set up.');
                return
            end
            fd = pmodel.FlowDiscretization;
            fp = pmodel.FlowPropertyFunctions;
            % Replace existing properties with total flux variants
            fd = fd.setStateFunction('PhaseFlux', PhaseFluxFixedTotalVelocity(pmodel));
            fd = fd.setStateFunction('PhaseUpwindFlag', PhasePotentialUpwindFlag(pmodel));
            fd = fd.setStateFunction('ComponentPhaseFlux', ComponentPhaseFluxFractionalFlow(pmodel));
            % Set extra props
            fd = fd.setStateFunction('PhaseInterfacePressureDifferences', PhaseInterfacePressureDifferences(pmodel));
            fd = fd.setStateFunction('TotalFlux', FixedTotalFlux(pmodel));
            fd = fd.setStateFunction('FaceTotalMobility', FaceTotalMobility(pmodel));
            % Set flow properties
            if isTotalSat
                fp = fp.setStateFunction('TotalSaturation', TotalSaturation(pmodel));
                fp = fp.setStateFunction('ComponentPhaseDensity', ComponentPhaseDensityTotalSaturation(pmodel));
                fp = fp.setStateFunction('ComponentMobility', ComponentMobilityTotalSaturation(pmodel));
                fp = fp.setStateFunction('ComponentPhaseMass', ComponentPhaseMassTotalSaturation(pmodel));
            end
            % Replace object
            model.parentModel.FlowDiscretization = fd;
            model.parentModel.FlowPropertyFunctions = fp;
            if hasFacility
                % Disable primary variables in transport!
                model.parentModel.FacilityModel.primaryVariableSet = 'none';
                fdp = model.parentModel.FacilityModel.FacilityFlowDiscretization;
                qf = WellPhaseFluxTotalFixed(model.parentModel);
                fdp = fdp.setStateFunction('PhaseFlux', qf);
                model.parentModel.FacilityModel.FacilityFlowDiscretization = fdp;
                model.parentModel.FacilityModel.doPostUpdate = true;
            end
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@PhysicalModel(model, state0, state, dt, drivingForces, varargin{:});
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            isS = strcmpi(problem.primaryVariables, 'sT');
            if any(isS)
                state = model.updateStateFromIncrement(state, dx, problem, 'sT', inf, model.dsMaxTotal);
                state = model.capProperty(state, 'sT', 1e-8);
                dx = dx(~isS);
            end
            problem.primaryVariables = problem.primaryVariables(~isS);
            
            [state, report] = model.parentModel.updateState(state, problem, dx, drivingForces);
        end
        
        function  [convergence, values, names] = checkConvergence(model, problem)
            [convergence, values, names] = model.parentModel.checkConvergence(problem);
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
            if strcmpi(name, 'sT')
                fn = 'sT';
                index = ':';
            else
                [fn, index] = model.parentModel.getVariableField(name, varargin{:});
            end
        end
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            % Prepare state and model (temporarily) before solving a time-step
            [model, state] = prepareTimestep@WrapperModel(model, state, state0, dt, drivingForces);
        end

        function [state, report] = updateAfterConvergence(model, varargin)
            [state, report] = updateAfterConvergence@WrapperModel(model, varargin{:});
        end
    end
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
