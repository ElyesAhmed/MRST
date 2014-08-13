classdef PhysicalModel
%Base class for physical models
%
% SYNOPSIS:
%   model = PhysicalModel(G)
%
% DESCRIPTION:
%   Base class for implementing physical models for use with automatic
%   differentiation. This class cannot be used directly.
%
%   A physical model instance contains the functions for getting residuals
%   and jacobians, making a single nonlinear step and verifying
%   convergence. It also contains the functions for updating the state
%   based on the increments found by the linear solver so that the values
%   are physically correct.
%
% REQUIRED PARAMETERS:
%
%   G     - Simulation grid.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   See class properties.
%
% RETURNS:
%   Class instance.
%
% SEE ALSO:
%   ThreePhaseBlackOilModel, TwoPhaseOilWaterModel, ReservoirModel

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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

    properties
        % Operators used for construction of systems
        operators
        % Inf norm tolerance for nonlinear iterations
        nonlinearTolerance
        % Grid
        G
        % Verbosity from model routines
        verbose
    end
    
    methods
        function model = PhysicalModel(G, varargin) %#ok
            model.nonlinearTolerance = 1e-6;
            model.verbose = mrstVerbose();
            
            model = merge_options(model, varargin{:});
            
            % Physical model
            model.G = G;
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin) %#ok
            % Get the equations governing the system
            error('Base class not meant for direct use')
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments
            for i = 1:numel(problem.primaryVariables);
                 p = problem.primaryVariables{i};
                 % Update the state
                 state = model.updateStateFromIncrement(state, dx, problem, p);
            end
            report = [];
        end
        
        function state = updateAfterConvergence(model, state0, state, drivingForces) %#ok
            % Update state based on non-linear increment after timestep has
            % converged. Defaults to doing nothing since not all models
            % require this.
        end
        
        function [convergence, values] = checkConvergence(model, problem, n)
            % Check and report convergence based on residual tolerances
            if nargin == 2
                n = inf;
            end
            
            values = norm(problem, n);
            convergence = all(values < model.nonlinearTolerance);
            
            if model.verbose
                for i = 1:numel(values)
                    fprintf('%s (%s): %2.2e\t', problem.equationNames{i}, problem.types{i}, values(i));
                end
                fprintf('\n')
            end
        end
        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolve, nonlinsolve, onlyCheckConvergence, varargin)
            % Make a single linearized timestep
            [problem, state] = model.getEquations(state0, state, dt, drivingForces, ...
                                       'ResOnly', onlyCheckConvergence, varargin{:});
            [convergence, values] = model.checkConvergence(problem);
            
            % Defaults
            failureMsg = '';
            failure = false;
            [linearReport, updateReport] = deal(struct());
            if ~(convergence && ~onlyCheckConvergence)
                % Get increments for Newton solver
                [dx, ~, linearReport] = linsolve.solveLinearProblem(problem, model);
                
                % Let the non-linear solver decide what to do with the
                % increments to get the best convergence
                dx = nonlinsolve.stabilizeNewtonIncrements(problem, dx);
                
                % Finally update the state. The physical model knows which
                % properties are actually physically reasonable.
                [state, updateReport] = model.updateState(state, problem, dx, drivingForces);
                if any(cellfun(@(d) ~all(isfinite(d)), dx))
                    failure = true;
                    failureMsg = 'Linear solver produced non-finite values.';
                end
            end
            report = model.makeStepReport(...
                            'LinearSolver', linearReport, ...
                            'UpdateState',  updateReport, ...
                            'Failure',      failure, ...
                            'FailureMsg',   failureMsg, ...
                            'Converged',    convergence, ...
                            'Residuals',    values);
        end
        
        function report = makeStepReport(model, varargin) %#ok
            report = struct('LinearSolver', [], ...
                            'UpdateState',  [], ...
                            'Failure',      false, ...
                            'FailureMsg',   '', ...
                            'Converged',    false, ...
                            'Residuals',    []);
            report = merge_options(report, varargin{:});
        end
        
        function [gradient, result, report] = solveAdjoint(model, solver, getState,...
                                    getObjective, schedule, gradient, itNo, scaling)
            % Solve adjoints
            if nargin == 7
               scaling = struct('rate', 1, 'pressure', 1);
            end
            
            dt_steps = schedule.step.val;
            
            current = getState(itNo);
            before    = getState(itNo - 1);
            dt = dt_steps(itNo);
            
            lookupCtrl = @(step) schedule.control(schedule.step.control(step));
            [~, forces] = model.getDrivingForces(lookupCtrl(itNo));
            problem = model.getEquations(before, current, dt, forces, 'iteration', inf, 'scaling', scaling);
            
            if itNo < numel(dt_steps)
                after    = getState(itNo + 1);
                dt_next = dt_steps(itNo + 1);
                
                [~, forces_p] = model.getDrivingForces(lookupCtrl(itNo + 1));
                problem_p = model.getEquations(current, after, dt_next, forces_p,...
                                    'iteration', inf, 'reverseMode', true, 'scaling', scaling);
            else
                problem_p = [];
            end
            [gradient, result, rep] = solver.solveAdjointProblem(problem_p,...
                                        problem, gradient, getObjective(itNo), model);
            report = struct();
            report.Types = problem.types;
            report.LinearSolverReport = rep;
        end
        
        function [vararg, driving] = getDrivingForces(model, control) %#ok
            % Setup and pass on driving forces. Dummy version for base
            % class.
            vararg = {};
            driving = struct();
        end
        
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state). This
            % always result in an error, as this model knows of no variables.
            [fn, index] = deal([]);
            
            if isempty(index)
                error('PhysicalModel:UnknownVariable', ...
                    ['State variable ''', name, ''' is not known to this model']);
            end
        end
        
        function p = getProp(model, state, name)
            % Get a property based on the name
            [fn, index] = model.getVariableField(name);
            p = state.(fn)(:, index);
        end
        
        function varargout = getProps(model, state, varargin)
            % Get multiple properties based on the name
            varargout = cellfun(@(x) model.getProp(state, x), ...
                                varargin, 'UniformOutput', false);
        end
        
        function state = incrementProp(model, state, name, increment)
            % Increment property based on name
            [fn, index] = model.getVariableField(name);
            p = state.(fn)(:, index)  + increment;
            state.(fn)(:, index) = p;
        end
        
        function state = setProp(model, state, name, value)
            % Set property to given value based on name
            [fn, index] = model.getVariableField(name);
            state.(fn)(:, index) = value;
        end
        
        function dv = getIncrement(model, dx, problem, name)
            % Find increment in linearized problem with given name, or
            % output zero if not found
            isVar = problem.indexOfPrimaryVariable(name);
            if any(isVar)
                dv = dx{isVar};
            else
                dv = 0;
            end
        end
        
        function [state, val, val0] = updateStateFromIncrement(model, state, dx, problem, name, relchangemax, abschangemax)
            % Update a state, with optionally a maximum relative change
            % applied.
            if iscell(dx)
                dv = model.getIncrement(dx, problem, name);
            else
                % Numerical value, increment directly and do not safety
                % check that this is a part of the model
                dv = dx;
            end
            
            val0 = model.getProp(state, name);
            
            [changeRel, changeAbs] = deal(1);
            if nargin > 5
                [~, changeRel] = model.limitUpdateRelative(dv, val0, relchangemax);
            end
            if nargin > 6
                [~, changeAbs] = model.limitUpdateAbsolute(dv, abschangemax);
            end            
            % Limit update by lowest of the relative and absolute limits 
            change = min(changeAbs, changeRel);
            
            val   = val0 + dv.*repmat(change, 1, size(dv, 2));
            state = model.setProp(state, name, val);
        end
    end

    methods (Static)
        function [dv, change] = limitUpdateRelative(dv, val, maxRelCh)
            % Limit a update by relative limit
            biggestChange = max(abs(dv./val), [], 2);
            change = min(maxRelCh./biggestChange, 1);
            dv = dv.*repmat(change, 1, size(dv, 2));
        end
        
        function [dv, change] = limitUpdateAbsolute(dv, maxAbsCh)
            % Limit a update by absolute limit
            biggestChange = max(abs(dv), [], 2);
            change = min(maxAbsCh./biggestChange, 1);
            dv = dv.*repmat(change, 1, size(dv, 2));
        end
        
        function [vars, isRemoved] = stripVars(vars, names)
            isRemoved = cellfun(@(x) any(strcmpi(names, x)), vars);
            vars(isRemoved) = [];
        end
    end

end

