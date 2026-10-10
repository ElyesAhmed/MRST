classdef NashOptimizationTwinProblem < handle
    % NashOptimizationTwinProblem - Specialized for twin-problem Nash optimization
    % of coarse-grid models with parameter exchange
    
    properties 
        twinParameters            % Cell array of parameter configurations for each twin
        twinProblems              % Cell array of simulation problems for each twin
        twinReferenceStates       % Cell array of reference states for each twin
        twinObjectives            % Cell array of objective functions for each twin
        
        % Nash optimization parameters
        nashParameters = struct('maxIterations', 15, ...
            'convergenceTol', 1e-4, ...
            'mixingParameter', 0.5, ...
            'maxInnerIterations', 25, ...
            'verbose', true, ...
            'regularizationParam', [0.0,0.0], ...
            'priorWeights', []);
    end
    
    properties (SetAccess = protected)
        currentNashIteration = 0
        parameterMappings       % Mappings for parameter exchange
        optimizationHistory     % History of optimization
        optimalParams1          % Optimal parameters for twin 1
        optimalParams2          % Optimal parameters for twin 2
    end
    
    methods
        function obj = NashOptimizationTwinProblem(varargin)
            % Constructor for twin-problem Nash optimization
            
            % Use inputParser for robust argument handling
            ip = inputParser;

            
            % TwinProblem specific parameters
            ip.addParameter('twinParameters', {});
            ip.addParameter('twinProblems', {});
            ip.addParameter('twinReferenceStates', {});
            ip.addParameter('twinObjectives', {});
            ip.addParameter('nashParameters', struct());
            
            ip.parse(varargin{:});
            opt = ip.Results;
            
            % Store twin-specific properties
            obj.twinParameters = opt.twinParameters;
            obj.twinProblems = opt.twinProblems;
            obj.twinReferenceStates = opt.twinReferenceStates;
            obj.twinObjectives = opt.twinObjectives;
                    
            opt.nashParameters     = merge_options(opt.nashParameters, varargin{:});

            % Handle nashParameters
            if ~isempty(opt.nashParameters) && isstruct(opt.nashParameters)
                userParams = opt.nashParameters;
                if isfield(userParams, 'maxIterations')
                    obj.nashParameters.maxIterations = userParams.maxIterations;
                end
                if isfield(userParams, 'convergenceTol')
                    obj.nashParameters.convergenceTol = userParams.convergenceTol;
                end
                if isfield(userParams, 'mixingParameter')
                    obj.nashParameters.mixingParameter = userParams.mixingParameter;
                end
                if isfield(userParams, 'maxInnerIterations')
                    obj.nashParameters.maxInnerIterations = userParams.maxInnerIterations;
                end
                   
                if isfield(userParams, 'regularizationParam')
                    obj.nashParameters.regularizationParam = userParams.regularizationParam;
                end
                if isfield(userParams, 'verbose')
                    obj.nashParameters.verbose = userParams.verbose;
                end
            end
            
            % Initialize parameter mappings
            obj = obj.initializeParameterMappings();
        end
        
        function obj = initializeParameterMappings(obj)
            % Initialize mappings for parameter exchange between twins
            obj.parameterMappings = struct();
            
            if numel(obj.twinParameters) >= 2
                % Get parameter names for both twins - handle different sizes
                names1 = cell(1, numel(obj.twinParameters{1}));
                for i = 1:numel(obj.twinParameters{1})
                    names1{i} = obj.twinParameters{1}{i}.name;
                end
                
                names2 = cell(1, numel(obj.twinParameters{2}));
                for i = 1:numel(obj.twinParameters{2})
                    names2{i} = obj.twinParameters{2}{i}.name;
                end
                
                % Create mapping from twin1 to twin2
                obj.parameterMappings.twin1_to_twin2 = containers.Map();
                for i = 1:numel(names1)
                    idx = find(strcmp(names2, names1{i}));
                    if ~isempty(idx)
                        obj.parameterMappings.twin1_to_twin2(names1{i}) = idx;
                    end
                end
                
                % Create mapping from twin2 to twin1
                obj.parameterMappings.twin2_to_twin1 = containers.Map();
                for i = 1:numel(names2)
                    idx = find(strcmp(names1, names2{i}));
                    if ~isempty(idx)
                        obj.parameterMappings.twin2_to_twin1(names2{i}) = idx;
                    end
                end
            end
        end


        % Helper function to add regularization properly
        function [f, grad] = addRegularizationToObjective(obj, p, base_obj_fn, reg_weight, p0)
            % Call base objective to get [f_base, grad_base]
            if nargout > 1
                [f_base, grad_base] = base_obj_fn(p);
            else
                f_base = base_obj_fn(p);
                grad_base = [];
            end

            % Add regularization term
            reg_term = reg_weight * norm(p - p0)^2;
            f = f_base + reg_term;

            % Add regularization gradient if available
            if nargout > 1 && ~isempty(grad_base)
                grad_reg = 2 * reg_weight * (p - p0);
                grad = grad_base + grad_reg';
            else
                grad = [];
            end
        end

function [u_opt, history, nashResults] = optimizeNashEquilibrium(obj, initialGuess, varargin)
    % Optimize using Nash equilibrium for twin problems
    
    % Use inputParser for robust option handling
    ip = inputParser;
    ip.addParameter('maxIterations', obj.nashParameters.maxIterations);
    ip.addParameter('convergenceTol', obj.nashParameters.convergenceTol);
    ip.addParameter('mixingParameter', obj.nashParameters.mixingParameter);
    ip.addParameter('maxInnerIterations', obj.nashParameters.maxInnerIterations);
    ip.addParameter('verbose', obj.nashParameters.verbose);
    ip.addParameter('optimizer_1', 'unitboxbfgs');
    ip.addParameter('optimizer_2', 'unitboxlm');
    ip.addParameter('objectiveFunctions', {});
    
    ip.parse(varargin{:});
    opt = ip.Results;
    
    if nargin < 2 || isempty(initialGuess)
        error('Initial guess required for twin problem optimization');
    end
    
    % Initialize optimization variables
    nTwins = numel(obj.twinProblems);
    currentParams = cell(nTwins, 1);
    history = struct();
    nashResults = struct('iteration', [], 'objectiveValues', [], 'parameters', []);
    
    obj.currentNashIteration = 0;
    
    % Convert initial guess to appropriate format
    if iscell(initialGuess) && numel(initialGuess) == nTwins
        currentParams = initialGuess;
    else
        % Assume initialGuess is a combined vector, split it
        nParams = cellfun(@(p) sum(cellfun(@(x) x.nParam, p)), obj.twinParameters);
        currentParams{1} = initialGuess(1:nParams(1));
        if nTwins > 1
            currentParams{2} = initialGuess(nParams(1)+1:end);
        end
    end
    
  % In optimizeNashEquilibrium method:
if obj.nashParameters.regularizationParam(1) > 0
    for twinIdx = 1:nTwins
        reg_weight = obj.nashParameters.regularizationParam(twinIdx);
        obj_fns{twinIdx} = @(p) obj.addRegularizationToObjective(p, opt.objectiveFunctions{twinIdx}, reg_weight, currentParams{twinIdx});
    end
else
    obj_fns = opt.objectiveFunctions;
end
    
    % Main Nash optimization loop
    for nashIter = 1:opt.maxIterations
        obj.currentNashIteration = nashIter;
        
        if opt.verbose
            fprintf('\n=== Nash Twin-Problem Iteration %d/%d ===\n', ...
                nashIter, opt.maxIterations);
        end
        
        % Store previous parameters for convergence check
        previousParams = currentParams;
        
        % Optimize each twin sequentially
        for twinIdx = 1:nTwins
            if opt.verbose
                fprintf('Optimizing Twin %d/%d...\n', twinIdx, nTwins);
            end
            
            % Use the pre-defined objective function with unitBoxLM
            % Pass options as key/value pairs, not as a structure
            if strcmp(opt.optimizer_1, 'unitboxbfgs')

                [~, currentParams{twinIdx}, hist] = unitBoxBFGS(currentParams{twinIdx}, obj_fns{twinIdx},'objChangeTol', 1e-5, ...
                    'maxIt', opt.maxInnerIterations,  'logPlot', true);
            elseif strcmp(opt.optimizer_2, 'unitboxlm')

                                [~, currentParams{twinIdx}, hist] = unitBoxLM(currentParams{twinIdx}, obj_fns{twinIdx}, ...
                    'maxIt', opt.maxInnerIterations, 'resChangeTolRel', 1.0e-5);
            end
            
            % Store history
            if isempty(history)
                history = struct('twin', cell(1, nTwins));
            end
            history(nashIter).twin{twinIdx} = hist;
        end
        
        % Exchange parameters between twins
        if nTwins > 1
            currentParams = obj.exchangeParameters(currentParams, opt.mixingParameter);
        end
        
        % Store results
        objValues = arrayfun(@(twin) history(nashIter).twin{twin}.val(end), 1:nTwins);
        nashResults.iteration(nashIter) = nashIter;
        nashResults.objectiveValues(nashIter, :) = objValues;
        nashResults.parameters{nashIter} = currentParams;
        
        if opt.verbose
            fprintf('Objective values after iteration %d: ', nashIter);
            fprintf('%.4e ', objValues);
            fprintf('\n');
        end
        
        % Check convergence
        if nashIter > 1
            paramChanges = cellfun(@(curr, prev) norm(curr - prev) / (norm(prev) + eps), ...
                currentParams, previousParams);
            
            if opt.verbose
                fprintf('Parameter changes: ');
                fprintf('%.4e ', paramChanges);
                fprintf('\n');
            end
            
            if all(paramChanges < opt.convergenceTol)
                if opt.verbose
                    fprintf('Nash equilibrium converged after %d iterations\n', nashIter);
                end
                break;
            end
        end
    end
    
    % Combine optimal parameters
    u_opt = currentParams;
    
    % Store final results
    obj.optimizationHistory = history;
    obj.optimalParams1 = currentParams{1};
    if nTwins > 1
        obj.optimalParams2 = currentParams{2};
    end
end        
        function exchangedParams = exchangeParameters(obj, currentParams, mixingParameter)
            % Exchange parameters between twins with mixing
            
            exchangedParams = currentParams;
            
            if numel(currentParams) >= 2
                % Convert parameters to cell arrays for easy manipulation
                params1_cell = obj.parametersToCell(currentParams{1}, obj.twinParameters{1});
                params2_cell = obj.parametersToCell(currentParams{2}, obj.twinParameters{2});
                
                % Exchange parameters from twin1 to twin2
                if isfield(obj.parameterMappings, 'twin1_to_twin2')
                    keys = obj.parameterMappings.twin1_to_twin2.keys;
                    for i = 1:numel(keys)
                        paramName = keys{i};
                        twin2_idx = obj.parameterMappings.twin1_to_twin2(paramName);
                        
                        % Find corresponding parameter in twin1
                        twin1_idx = find(cellfun(@(x) strcmp(x.name, paramName), obj.twinParameters{1}));
                        
                        if ~isempty(twin1_idx) && twin2_idx <= numel(obj.twinParameters{2})
                            % Unscale and mix parameters
                            pu1 = obj.twinParameters{1}{twin1_idx}.unscale(params1_cell{twin1_idx});
                            pu2 = obj.twinParameters{2}{twin2_idx}.unscale(params2_cell{twin2_idx});
                            
                            pu_mixed = mixingParameter * pu1 + (1 - mixingParameter) * pu2;
                            
                            % Scale back and update
                            params2_cell{twin2_idx} = obj.twinParameters{2}{twin2_idx}.scale(pu_mixed);
                        end
                    end
                end
                
                % Exchange parameters from twin2 to twin1
                if isfield(obj.parameterMappings, 'twin2_to_twin1')
                    keys = obj.parameterMappings.twin2_to_twin1.keys;
                    for i = 1:numel(keys)
                        paramName = keys{i};
                        twin1_idx = obj.parameterMappings.twin2_to_twin1(paramName);
                        
                        % Find corresponding parameter in twin2
                        twin2_idx = find(cellfun(@(x) strcmp(x.name, paramName), obj.twinParameters{2}));
                        
                        if ~isempty(twin2_idx) && twin1_idx <= numel(obj.twinParameters{1})
                            % Unscale and mix parameters
                            pu2 = obj.twinParameters{2}{twin2_idx}.unscale(params2_cell{twin2_idx});
                            pu1 = obj.twinParameters{1}{twin1_idx}.unscale(params1_cell{twin1_idx});
                            
                            pu_mixed = mixingParameter * pu2 + (1 - mixingParameter) * pu1;
                            
                            % Scale back and update
                            params1_cell{twin1_idx} = obj.twinParameters{1}{twin1_idx}.scale(pu_mixed);
                        end
                    end
                end
                
                % Convert back to vectors
                exchangedParams{1} = vertcat(params1_cell{:});
                exchangedParams{2} = vertcat(params2_cell{:});
            end
        end
        
        function param_cell = parametersToCell(~, param_vector, parameters)
            % Convert parameter vector to cell array based on parameter structure
            nparam = cellfun(@(x) x.nParam, parameters);
            param_cell = mat2cell(param_vector, nparam, 1);
        end
        
        function plotTwinConvergence(obj)
            % Plot convergence history for twin problems
            if isempty(obj.optimizationHistory)
                error('No optimization history available. Run optimizeNashEquilibrium() first.');
            end
            
            nashIters = numel(obj.optimizationHistory);
            objValues = zeros(nashIters, numel(obj.twinProblems));
            
            for i = 1:nashIters
                for j = 1:numel(obj.twinProblems)
                    objValues(i, j) = obj.optimizationHistory(i).twin{j}.val(end);
                end
            end
            
            figure;
            subplot(2,1,1);
            plot(1:nashIters, objValues, 'LineWidth', 2);
            xlabel('Nash Iteration');
            ylabel('Objective Value');
            title('Twin-Problem Nash Equilibrium Convergence');
            legend(arrayfun(@(x) sprintf('Twin %d', x), 1:numel(obj.twinProblems), 'UniformOutput', false));
            grid on;
            
            subplot(2,1,2);
            semilogy(1:nashIters, abs(diff([zeros(1, size(objValues, 2)); objValues])), '--', 'LineWidth', 2);
            xlabel('Nash Iteration');
            ylabel('Absolute Improvement');
            title('Objective Improvement per Iteration');
            legend(arrayfun(@(x) sprintf('Twin %d', x), 1:numel(obj.twinProblems), 'UniformOutput', false));
            grid on;
        end
        
        function [setup1, setup2] = getOptimalSetups(obj)
            % Get optimal setups after optimization
            if isempty(obj.optimizationHistory)
                error('Run optimizeNashEquilibrium() first to get optimal parameters');
            end
            
            setup1 = updateSetupFromScaledParameters(...
                obj.twinProblems{1}.SimulatorSetup, obj.twinParameters{1}, obj.optimalParams1);
            
            if numel(obj.twinProblems) > 1
                setup2 = updateSetupFromScaledParameters(...
                    obj.twinProblems{2}.SimulatorSetup, obj.twinParameters{2}, obj.optimalParams2);
            else
                setup2 = [];
            end
        end
    end
end