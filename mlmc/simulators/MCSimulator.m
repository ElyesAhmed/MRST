classdef MCSimulator < MRSTEnsemble
    
    properties(Access = protected)
        estimate   = 0
        variance   = 0
        cost       = 0
        numSamples = 0
        rmse       = inf
        history    = {}
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function mc = MCSimulator(mrstExample, samples, qoi, varargin)
            assert(numel(qoi.names) == 1, 'MCSimulator only supports one qoi');
            mc = mc@MRSTEnsemble(mrstExample, samples, qoi, varargin{:});
            mc.figures.progressMC = [];
        end
        
        %-----------------------------------------------------------------%
        function runSimulation(mc, varargin)
            opt = struct('range'       , inf , ...
                         'batchSize'   , inf , ...
                         'maxSamples'  , 10  , ...
                         'minSamples'  , 2   , ...
                         'tolerance'   , -inf, ...
                         'relTolerance', 1e-2, ...
                         'prompt'      , true, ...
                         'plotProgress', true);
            opt = merge_options(opt, varargin{:});
            assert(xor(isfinite(opt.tolerance), isfinite(opt.relTolerance)), ...
                   'Choose either relative or absolute tolerance');
            % Check for computed samples
            range = mc.getComputedSampleRange();
            n     = numel(range);
            if n > 0 && opt.prompt
                % We found computed samples - ask user if they should be
                % included in the estimate
                prompt = sprintf(['Found %d QoIs for this problem. '   , ...
                                  'Would you like to include these in ', ...
                                  'your estimate (y)? If not, I will ' , ...
                                  'delete them y/n [y]: '], n          );
                str = input(prompt,'s');
                if strcmpi(str, 'y') || isempty(str)
                    % Include samples in estimate
                    fprintf('Ok, I will include the QoIs in the estimate.\n');
                    if mc.numSamples ~= n
                        mc.resetStatistics();
                        mc.updateStatistics(range);
                    end
                else
                    % Reset simulator. This includes wiping out any samples
                    % already computed
                    mc.reset('prompt', false);
                end
            end
            % Set tolerance
            tolerance = opt.tolerance;
            if isinf(opt.tolerance)
                % Set tolerance based on relative reduction target
                tolerance = opt.relTolerance*mc.estimate;
            end
            % Compute number of samples needed
            n = mc.computeNumSamples(tolerance, opt);
            while (any(n > 0) && mc.numSamples < opt.maxSamples) ... 
                              || mc.numSamples < opt.minSamples
                % Run batch of n new samples
                range = mc.runBatch('batchSize', n, 'plotProgress', opt.plotProgress);
                % Update statistics
                mc.updateStatistics(range);
                if isinf(opt.tolerance)
                    % Adjust tolerance if we aim at relative reduction
                    tolerance = opt.relTolerance*mc.estimate;
                end
                % Recompute number of samples needed
                n = mc.computeNumSamples(tolerance, opt);
                if opt.plotProgress
                    % Plot progress of estimate with rmse bounds
                    out = mc.getHistory();
                    mc.figures.progressMC = plotMonteCarloProgress(out, ...
                                                    mc.figures.progressMC);
                    drawnow(); pause(0.1);
                end
            end
            
        end
        
        %-----------------------------------------------------------------%
        function range = runBatch(mc, varargin)
            range = mc.simulateEnsembleMembers(varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function maxId = getLargestSeed(mc)
            maxId = max(mc.qoi.ResultHandler.getValidIds());
            if isempty(maxId), maxId = 0; end
        end
        
        %-----------------------------------------------------------------%
        function range = getComputedSampleRange(mc)
            range = mc.qoi.ResultHandler.getValidIds();
        end
        
        %-----------------------------------------------------------------%
        function reset(mc, varargin)
            reset@MRSTEnsemble(mc, varargin{:});
            mc.resetStatistics();
        end
        
        %-----------------------------------------------------------------%
        function resetStatistics(mc)
            mc.estimate   = 0;
            mc.variance   = 0;
            mc.cost       = 0;
            mc.numSamples = 0;
            mc.rmse       = inf;
            mc.history    = {};
        end
        
        %-----------------------------------------------------------------%
        function updateStatistics(mc, range)
            if nargin < 2 || all(isinf(range))
                range = mc.qoi.ResultHandler.getValidIds();
                mc.resetStatistics();
            end
            % Get current statistics
            m0 = mc.estimate;
            v0 = mc.variance;
            c0 = mc.cost;
            n0 = mc.numSamples;
            % Get batch statistics
            [m, v] = mc.qoi.computeQoIMean(range);
            c = m.cost; m = m.(mc.qoi.names{1}); v = v.(mc.qoi.names{1});
            n = numel(range);
            % Update statistics
            mc.estimate   = computeMean(m0, m, n0, n);
            mc.variance   = computeVariance(v0, v, m0, m, n0, n);
            mc.cost       = computeMean(c0, c, n0, n);
            mc.numSamples = mc.numSamples + numel(range);
            if mc.numSamples > 1
                mc.rmse = sqrt(mc.variance/mc.numSamples);
            end
            % Update history
            mc.history{end+1} = struct('estimate'  , mc.estimate  , ...
                                       'variance'  , mc.variance  , ...
                                       'cost'      , mc.cost      , ...
                                       'numSamples', mc.numSamples, ...
                                       'rmse'      , mc.rmse      );
        end
        
        %-----------------------------------------------------------------%
        function n = computeNumSamples(mc, tolerance, opt)
            n = ceil(mc.variance/tolerance^2) - mc.numSamples;
            n = min(n, opt.maxSamples - mc.numSamples);
            n = max(n, opt.minSamples - mc.numSamples);
            n = min(n, opt.batchSize);
            n = max(n, 0);
        end
        
        %-----------------------------------------------------------------%
        % Getters                                                         %
        %-----------------------------------------------------------------%
        function stat = getStatistics(mc)
            stat = struct();
            stat.estimate   = mc.estimate;
            stat.variance   = mc.variance;
            stat.cost       = mc.cost;
            stat.numSamples = mc.numSamples;
            stat.rmse       = mc.rmse;
        end
        
        %-----------------------------------------------------------------%
        function history = getHistory(mc)
            get = @(fn) reshape(cellfun(@(h) h.(fn), mc.history), [], 1);
            history = struct();
            history.estimate   = get('estimate');
            history.variance   = get('variance');
            history.cost       = get('cost');
            history.numSamples = get('numSamples');
            history.rmse       = get('rmse');
        end

    end
    
end