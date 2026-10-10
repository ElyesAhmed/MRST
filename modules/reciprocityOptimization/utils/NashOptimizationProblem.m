classdef NashOptimizationProblem < OptimizationProblem
    properties 
        objectives                  % cell array of objective functions or structs
        controlTypes               % cell array of control types for each objective
        nashParameters             % parameters for Nash equilibrium optimization
    end
    properties (SetAccess = protected)
        fullMaps                    % Store the original full mappings
        fullUInit                   % Store the original full initial controls
        fullBounds                  % Store the original full bounds
        objectiveMaps               % cell array of mappings for each objective
        objectiveBounds             % cell array of bounds for each objective
        currentNashIteration = 0    % Track current Nash iteration
    end
    
    methods
        function p = NashOptimizationProblem(samples, varargin)
            % Use inputParser for robust argument handling
            ip = inputParser;
            
            % BaseEnsemble parameters
            ip.addParameter('directory', '');
            ip.addParameter('name', 'tmp');
            
            % OptimizationProblem parameters
            ip.addParameter('objective', []);
            ip.addParameter('constraints', []);
            ip.addParameter('solverFun', []);
            ip.addParameter('maps', []);
            ip.addParameter('parameters', []);
            ip.addParameter('bounds', []);
            ip.addParameter('statFun', @mean);
            ip.addParameter('updateProblemFun', []);
            ip.addParameter('fetchVariablesFun', []);
            ip.addParameter('plotProgress', true);
            ip.addParameter('setupType', '');
            ip.addParameter('background', false);
            ip.addParameter('regularization', []);
            ip.addParameter('clearStates', false);
            ip.addParameter('solverFunOptions', {});
            
            % NashOptimizationProblem specific parameters
            ip.addParameter('objectives', []);
            ip.addParameter('controlTypes', {});
            ip.addParameter('nashParameters', struct());
            ip.addParameter('u_init', []);
            
            ip.parse(varargin{:});
            opt = ip.Results;
            
            % Handle directory
            if isempty(opt.directory)
                opt.directory = fullfile(mrstOutputDirectory(), opt.name);
            end
            
            % Call parent constructor
            parentArgs = {};
            parentFields = {'directory', 'name', 'objective', 'constraints', ...
                'solverFun', 'maps', 'parameters', 'bounds', 'statFun', ...
                'updateProblemFun', 'fetchVariablesFun', 'plotProgress', ...
                'setupType', 'regularization', 'clearStates', 'solverFunOptions'};
            
            for i = 1:numel(parentFields)
                if isfield(opt, parentFields{i})
                    parentArgs = [parentArgs, parentFields{i}, {opt.(parentFields{i})}];
                end
            end
            
            p = p@OptimizationProblem(samples, parentArgs{:});
            
            % Store Nash-specific properties
            p.objectives = opt.objectives;
            p.controlTypes = opt.controlTypes;
            p.fullMaps = opt.maps;
            p.fullUInit = opt.u_init;
            p.fullBounds = opt.maps.bounds;
            
            % Default Nash parameters
            defaultNashParams = struct('maxIterations', 10, ...
                'convergenceTol', 1e-4, ...
                'mixingParameter', 0.5, ...
                'verbose', true);
            
            % Handle nashParameters
            p.nashParameters = defaultNashParams;
            if ~isempty(opt.nashParameters) && isstruct(opt.nashParameters)
                userParams = opt.nashParameters;
                if isfield(userParams, 'maxIterations')
                    p.nashParameters.maxIterations = userParams.maxIterations;
                end
                if isfield(userParams, 'convergenceTol')
                    p.nashParameters.convergenceTol = userParams.convergenceTol;
                end
                if isfield(userParams, 'mixingParameter')
                    p.nashParameters.mixingParameter = userParams.mixingParameter;
                end
                if isfield(userParams, 'verbose')
                    p.nashParameters.verbose = userParams.verbose;
                end
            end
            
            % Create reduced mappings for each objective
            p = p.createObjectiveMappings();
        end
        
        function p = createObjectiveMappings(p)
            % Create reduced mappings and bounds for each objective
            nObjectives = numel(p.objectives);
            p.objectiveMaps = cell(nObjectives, 1);
            p.objectiveBounds = cell(nObjectives, 1);
            
            for objIdx = 1:nObjectives
                if ~isempty(p.controlTypes) && numel(p.controlTypes) >= objIdx && ~isempty(p.controlTypes{objIdx})
                    % Use reduceControlMappings for specific control types
                    [p.objectiveMaps{objIdx}, ~] = reduceControlMappings(...
                        p.fullMaps, p.fullUInit, 'type', p.controlTypes{objIdx});
                    
                    % Create corresponding bounds
                    if ~isempty(p.fullBounds)
                        fullControlIndices = [];
                        for j = 1:numel(p.fullMaps.type)
                            for i = 1:numel(p.objectiveMaps{objIdx}.type)
                                if strcmp(p.objectiveMaps{objIdx}.type{i}, p.fullMaps.type{j}) && ...
                                        p.objectiveMaps{objIdx}.stepNo(i) == p.fullMaps.stepNo(j) && ...
                                        p.objectiveMaps{objIdx}.wellNo(i) == p.fullMaps.wellNo(j)
                                    fullControlIndices(end+1) = j;
                                    break;
                                end
                            end
                        end
                        p.objectiveBounds{objIdx} = p.fullBounds(fullControlIndices, :);
                    else
                        p.objectiveBounds{objIdx} = [];
                    end
                else
                    % Default: use all controls
                    p.objectiveMaps{objIdx} = p.fullMaps;
                    p.objectiveBounds{objIdx} = p.fullBounds;
                end
            end
        end
        
        function varargout = getScaledObjective(p, us, objIndex, isConstraint)
            if nargin < 4
                isConstraint = false;
            end
            if nargin < 3
                objIndex = 1;
            end
            
            % Use the appropriate mapping and bounds for this objective
            originalMaps = p.maps;
            originalBounds = p.bounds;
            p.maps = p.objectiveMaps{objIndex};
            p.bounds = p.objectiveBounds{objIndex};
            p.objective = p.objectives{objIndex};
            
            % Call parent method
            if nargout < 2
                varargout{1} = getScaledObjective@OptimizationProblem(p, us, isConstraint);
            else
                [varargout{1}, varargout{2}] = getScaledObjective@OptimizationProblem(p, us, isConstraint);
            end
            
            % Restore original mapping and bounds
            p.maps = originalMaps;
            p.bounds = originalBounds;
        end
        
        function [u_opt, history, nashResults] = optimizeNashEquilibrium(p, initialGuess, varargin)
            opt = struct('maxIterations', p.nashParameters.maxIterations, ...
                        'convergenceTol', p.nashParameters.convergenceTol, ...
                        'mixingParameter', p.nashParameters.mixingParameter, ...
                        'verbose', p.nashParameters.verbose, ...
                        'optimizerOptions', {{}});
            [opt, optimopts] = merge_options(opt, varargin{:});
            
            if nargin < 2 || isempty(initialGuess)
                ids = p.iterationControls.getValidIds();
                if ~isempty(ids) && ids(1)==1
                    us0 = p.iterationControls{1};
                else
                    error('Initial guess required');
                end
            else
                if ~isa(initialGuess, 'double')
                    try
                        initialGuess = p.fetchVariablesFun(initialGuess);
                    catch
                        error('Not able to fetch decision variables from initialGuess');
                    end
                end
                us0 = p.scaleVariables(initialGuess);
            end
            
            nObjectives = numel(p.objectives);
            currentControls = cell(nObjectives, 1);
            history = struct();
            nashResults = struct('iteration', [], 'objectiveValues', [], 'controls', []);
            
            p.currentNashIteration = 0;
            
            % Convert full controls to each objective's reduced space
            fullControls = us0;
            for objIdx = 1:nObjectives
                currentControls{objIdx} = p.reduceToObjectiveSpace(fullControls, objIdx);
            end
            
            for nashIter = 1:opt.maxIterations
                p.currentNashIteration = nashIter;
                if opt.verbose
                    fprintf('\n--- Nash Iteration %d/%d ---\n', nashIter, opt.maxIterations);
                end
                
                % Store current controls for all objectives
                previousControls = currentControls;
                
                % Optimize each objective sequentially in its own control space
                for objIdx = 1:nObjectives
                    if opt.verbose
                        fprintf('Optimizing objective %d/%d...\n', objIdx, nObjectives);
                    end

                    p.reset('prompt', false)
                    % Optimize this objective in its reduced control space
                                % Use this objective's specific bounds
                    p.bounds = p.objectiveBounds{objIdx};
                    p.maps = p.objectiveMaps{objIdx};
                    [u_opt_obj, hist_obj] = p.optimizeSingleObjective(...
                        currentControls{objIdx}, objIdx, optimopts{:});
                    p.updateProblemFromIteraion(u_opt_obj, hist_obj, opt.mixingParameter);

                    
                    % Store optimized controls for this objective
                    currentControls{objIdx} = opt.mixingParameter.*u_opt_obj+(1-opt.mixingParameter).*hist_obj.u{1};
                    history(nashIter).objective(objIdx) = hist_obj;
                end
                
                % Combine controls from all objectives back to full space
                
                % Evaluate all objectives at current combined controls
                objValues = zeros(1, nObjectives);
                for objIdx = 1:nObjectives
                    % Get controls in this objective's reduced space
%                     reducedControls = p.reduceToObjectiveSpace(fullControls, objIdx);
%                     objVal = p.getScaledObjective(reducedControls, objIdx);
                    objValues(objIdx) = history(nashIter).objective(objIdx).val(end);
                end
                fullControls = [history(nashIter).objective(1).u{end} history(nashIter).objective(2).u{end}];
                % Store results
                nashResults.iteration(nashIter) = nashIter;
                nashResults.objectiveValues(nashIter, :) = objValues;
                nashResults.controls{nashIter} = fullControls;
                
                if opt.verbose
                    fprintf('Objective values after iteration %d: ', nashIter);
                    fprintf('%.4f ', objValues);
                    fprintf('\n');
                end
                
                % Check convergence
                if nashIter > 1
                    improvement = norm(objValues - nashResults.objectiveValues(nashIter-1, :)) / ...
                                 norm(nashResults.objectiveValues(nashIter-1, :));
                    if opt.verbose
                        fprintf('Relative improvement: %.6f\n', improvement);
                    end
                    if improvement < opt.convergenceTol
                        if opt.verbose
                            fprintf('Nash equilibrium converged after %d iterations\n', nashIter);
                        end
                        break;
                    end
                end
                close all;
            end
            
            u_opt = p.unscaleVariables(fullControls(:,objIdx));
        end
        
        function [u_opt, history] = optimizeSingleObjective(p, initialGuess, objIndex, varargin)
            opt = struct('maximize', true, ...
                        'optimizer', 'unitboxbfgs', ...
                         'background',     false, ...
                        'regularizationParam', 0);
            [opt, optimopts] = merge_options(opt, varargin{:});
            
            if opt.regularizationParam > 0
                p.regularization = getRegularizationFun(initialGuess, opt.regularizationParam);
            else
                p.regularization = [];
            end

            
            switch lower(opt.optimizer)
                case {'unitboxbfgs'}
                    [~, u_opt, history] = unitBoxBFGS(initialGuess, ...
                        @(u)p.getScaledObjective(u, objIndex), ...
                        optimopts{:}, 'maximize', opt.maximize);
                    %unitBoxBFGS(us0, @p.getScaledObjective, optimopts{:}, 'maximize', opt.maximize);
                    
                case 'unitboxlm'
                    [~, u_opt, history] = unitBoxLM(initialGuess, ...
                        @(u)p.getScaledObjective(u, objIndex), ...
                        optimopts{:});
                    
                otherwise
                    error('Unsupported optimizer: %s', opt.optimizer);
            end
            
        end
        
        function reducedControls = reduceToObjectiveSpace(p, fullControls, objIndex)
            % Convert full controls to reduced space for specific objective
            objMaps = p.objectiveMaps{objIndex};
            
            % For each control in the reduced mapping, find corresponding control in full mapping
            reducedControls = zeros(numel(objMaps.type), 1);
            
            for i = 1:numel(objMaps.type)
                % Find the corresponding control in the full mapping
                for j = 1:numel(p.fullMaps.type)
                    if strcmp(objMaps.type{i}, p.fullMaps.type{j}) && ...
                       objMaps.stepNo(i) == p.fullMaps.stepNo(j) && ...
                       objMaps.wellNo(i) == p.fullMaps.wellNo(j)
                        reducedControls(i) = fullControls(j);
                        break;
                    end
                end
            end
        end
        
        function updateProblemFromIteraion(p, uopt, history_opt, mixingParameter,varargin)

            ns =  p.samples.num;
            opt = struct('sample', []);
            opt = merge_options(opt, varargin{:});
            for seed=1:ns
                problem_or{seed} = p.getSampleProblem(seed, 'sample', opt.sample);
                 uc =  p.unscaleVariables(history_opt.u{1});    
                 p.samples.problem{seed} = p.updateProblemFun(problem_or{seed}, mixingParameter*p.unscaleVariables(uopt)+(1-mixingParameter)*uc);  
            end
        end
        
        function plotNashConvergence(p, nashResults)
            % Plot Nash equilibrium convergence
            figure;
            subplot(2,1,1);
            plot(nashResults.iteration, nashResults.objectiveValues, 'LineWidth', 2);
            xlabel('Nash Iteration');
            ylabel('Objective Value');
            title('Nash Equilibrium Convergence');
            legend(arrayfun(@(x)sprintf('Objective %d', x), 1:size(nashResults.objectiveValues, 2), 'UniformOutput', false));
            grid on;
            
            subplot(2,1,2);
            plot(nashResults.iteration, sum(nashResults.objectiveValues, 2), 'g-o', 'LineWidth', 2);
            xlabel('Nash Iteration');
            ylabel('Combined Objective Value');
            title('Total Objective Convergence');
            grid on;
        end
    end
end


% -----------------------------------------------------------------
% Helper functions for well-controls (non-empty 'maps'-prop)
% -----------------------------------------------------------------

function problem = updateWellControls(problem, u, maps, W)
schedule = problem.SimulatorSetup.schedule;
if isempty(maps)
    nc = 0;
else
    nc = numel(maps.type);
end
for k = 1:nc
    [sno, wno, tp] = deal(maps.stepNo(k), maps.wellNo(k), maps.type{k});
    if maps.isTarget(k)
        assert(strcmp(tp, schedule.control(sno).W(wno).type), ...
            'Schedule contains non-expected well type');
        schedule.control(sno).W(wno).val = u(k);
        % also update if part of lims
        if isfield(schedule.control(sno).W(wno).lims, tp)
            schedule.control(sno).W(wno).lims.(tp) = u(k);
        end
    else
        schedule.control(sno).W(wno).lims.(tp) = u(k);
    end
end
%problem.SimulatorSetup.schedule = schedule;
if isfield(W, 'posControl') && nc < numel(u)
    % trajectory controls
    % W = schedule.control(1).W;
    pcix = arrayfun(@(w)~isempty(w.posControl), W);
    pc = {W(pcix).posControl};
    ix = nc;
    % update interpolation points
    for k = 1:numel(pc)
        nParam = pc{k}.parameters.nParam;
        pc{k} =  pc{k}.control2param(u(ix + (1:nParam)));
        ix = ix + nParam;
    end
    gboModel =  problem.SimulatorSetup.model;
    if ~isa(gboModel, 'GenericBlackOilModel')
        gboModel = problem.SimulatorSetup.model.getReservoirModel; 
    end
    % update trajectories
    for ks = 1:numel(schedule.control)
        for kw = 1:numel(schedule.control(ks).W)
            if ~isempty(W(kw).posControl)
                schedule.control(ks).W(kw) = ...
                    updateWellTrajectory(gboModel, ...
                        schedule.control(ks).W(kw), [], ...
                        W(kw).posControl.getTrajectory);
                % make sure posControl is included in schedule
                schedule.control(ks).W(kw).posControl = W(kw).posControl;
            end
        end
    end
end
problem.SimulatorSetup.schedule = schedule;
end
% -----------------------------------------------------------------

function u = fetchWellControlValues(problem, maps, W)
if isempty(maps)
    u = [];
else
    u = zeros(numel(maps.type), 1);
end
schedule = problem.SimulatorSetup.schedule;
for k = 1:numel(u)
    [sno, wno, tp] = deal(maps.stepNo(k), maps.wellNo(k), maps.type{k});
    if maps.isTarget(k)
        assert(strcmp(tp, schedule.control(sno).W(wno).type), ...
            'Schedule contains non-expected well type');
        u(k) = schedule.control(sno).W(wno).val;
    else
        u(k) = schedule.control(sno).W(wno).lims.(tp);
    end
end

if isfield(W, 'posControl')
    % position is same over all steps
    %W = schedule.control(1).W;
    ix = arrayfun(@(w)~isempty(w.posControl), W);
    pc = {W(ix).posControl};
    up = cellfun(@(x)x.param2control(), pc, 'UniformOutput', false);
    u = [u; vertcat(up{:})];
end
end

% -------------------------------------------------------------------------
% Helper functions for parameter-problems (non-empty 'parameters'-prop)
% -------------------------------------------------------------------------

function problem = updateModelParameters(problem, u, params)
nparam = cellfun(@(x)x.nParam, params);
u = mat2cell(u, nparam, 1);
% Create new setup, and set parameter values
setup = problem.SimulatorSetup;
% make sure to recreate discretizations
setup.model = setup.model.removeStateFunctionGroupings();

for k = 1:numel(params)
    tmp   = params{k}.unscale(u{k});
    setup = params{k}.setParameter(setup, tmp);
end
problem.SimulatorSetup = setup;
end
% -------------------------------------------------------------------------

function u = fetchParameterValues(problem, params)
u = getScaledParameterVector(problem.SimulatorSetup, params);
end


% -------------------------------------------------------------------------
% Setup of default simulation/diagnsotics problems
% -------------------------------------------------------------------------
function p = setupProblem(p, opt)
props = intersect(propertynames(p), fieldnames(opt));
for k = 1:numel(props)
    p.(props{k}) = opt.(props{k});
end

if ~isempty(opt.setupType)
    % setup simulation/diagnostics problem
    isControlType     = ~isempty(p.maps);
    isParamType    = ~isempty(p.parameters);
    isPositionType = ~isempty(p.W);
    assert(isControlType || isParamType || isPositionType, ...
        'Default setup requires either input of ''maps'' or ''parameters''');
    if (isControlType || isPositionType) && isParamType
        warning('Can''t work with both well controls and parameters. %s', ...
                'Parameters will be disregarded');
        p.parameters = [];
    end
    if isempty(p.solverFun)
        if strcmp(opt.setupType, 'simulation')
            p.solverFun = @(problem, obj, h, computeGradient) ...
                simulationSolverFun(problem, obj, 'objectiveHandler', h, ...
                'computeGradient', computeGradient, ...
                'parameters', p.parameters, 'maps', p.maps, ...
                'clearStatesAfterAdjoint', opt.clearStates, opt.solverFunOptions{:});
        elseif strcmp(opt.setupType, 'diagnostics')
            p.solverFun = @(problem, obj, h, computeGradient) ...
                diagnosticsSolverFun(problem, obj, 'objectiveHandler', h, ...
                'computeGradient', computeGradient, ...
                'parameters', p.parameters, 'maps', p.maps, 'W', p.W, opt.solverFunOptions{:});
        else
            error('Unknown setup-type: %s', opt.setupType);
        end
    end
    if isempty(p.updateProblemFun)
        if isControlType || isPositionType   
            p.updateProblemFun = ...
                @(problem,u)updateWellControls(problem, u, p.maps, p.W);
        else
            p.updateProblemFun = ...
                @(problem,u)updateModelParameters(problem, u, p.parameters);
        end
    end
    if isempty(p.fetchVariablesFun)
        if isControlType || isPositionType
            p.fetchVariablesFun = ...
                @(problem)fetchWellControlValues(problem, p.maps, p.W);
        else% parameter problem
            p.fetchVariablesFun = ...
                @(problem)fetchParameterValues(problem, p.parameters);
        end
    end
    % deal with bounds
    p.bounds = zeros(2,0);
    if isControlType
        if ~isempty(p.bounds)
            warning('Disregarding input: ''bounds'' (using maps.bounds)');
        end
        p.bounds = p.maps.bounds;
    end
    if isPositionType
        if ~isempty(p.bounds)
            warning('Disregarding input: ''bounds'' (using pc.boxes)');
        end
        assert(isfield(p.W(1), 'posControl'));
        pc = vertcat(p.W.posControl);
        boxes = [pc.boxes];
        p.bounds = [p.bounds; horzcat(boxes{:})'];
    end
    if isParamType
        if ~isempty(p.bounds)
            warning('Disregarding input: ''bounds'' (using parameter.boxLims)');
        end
        p.bounds = [];
    end
end


end

% -------------------------------------------------------------------------
% Handler utilities
% -------------------------------------------------------------------------

function h = getHandler(maindd, prefix, controlNo, memberNo, showWarning)
if nargin < 5
    showWarning = true;
end
pth = fullfile(maindd, num2str(controlNo), num2str(memberNo));
if isfolder(pth)
    [dd, fld] = fileparts(pth);
    h = ResultHandler('dataPrefix', prefix, ...
        'writeToDisk', true,...
        'dataDirectory', dd, ...
        'dataFolder', fld, ...
        'cleardir', false);
else
    h = [];
    if showWarning
        warning('Non-existing directory: %s', pth)
    end
end
end

function problem = resetHandlers(problem, directory, folder)
if isfield(problem, 'OutputHandlers')
    fn = fieldnames(problem.OutputHandlers);
    for k = 1:numel(fn)
        h = problem.OutputHandlers.(fn{k});
        if isa(h, 'ResultHandler')
            if ~exist(fullfile(directory, folder), 'dir')
                mkdir(fullfile(directory, folder));
            end
            h.dataDirectory = directory;
            h.dataFolder    = folder;
            problem.OutputHandlers.(fn{k}) = h;
        end
    end
end
end
    
% -------------------------------------------------------------------------
% Regularization
% -------------------------------------------------------------------------

function fn = getRegularizationFun(us0, alpha)
fn = @(us)regularizationFn(us, us0, alpha);
end

function varargout = regularizationFn(us, us0, alpha)
varargout{1} = alpha*sum( (us-us0).^2 );
if nargout > 1
    varargout{2} = 2*alpha*(us-us0);
end
end


function x = applyStatFun(statfun, c)
% Apply some 'statistic' function to cell-array c, but also accept @vertcat
if strcmp(func2str(statfun), 'vertcat')
    x = statfun(c{:});
else
    [nRow, nCol] = size(c{1});
    if nCol > 1
        c = applyFunction(@(x)x(:), c);
    end
    x = statfun(horzcat(c{:}), 2);
    x = reshape(x, [nRow, nCol]);      
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
