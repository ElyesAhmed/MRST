function [misfitVals ,varargout] = evaluateMatchRegularizedSummands(pvec, obj, setup, parameters, states_ref, varargin)
% [Your existing header comments...]

opt = struct('Verbose',           mrstVerbose(),...
             'NonlinearSolver', [],...
             'objScaling',1, ...
             'enforceBounds',  true, ...
             'accumulateResiduals', [],...
             'regularizationParam', 0.0, ...        % ADD THIS
             'initialParameters', []);              % ADD THIS

[opt, extra] = merge_options(opt, varargin{:});

% ADD REGULARIZATION CALCULATION HERE
reg_term = 0;
reg_grad = [];
if opt.regularizationParam > 0 && ~isempty(opt.initialParameters)
    reg_weight = opt.regularizationParam;
    reg_term = reg_weight * norm(pvec - opt.initialParameters)^2;
    
    % Compute regularization gradient if needed
    if nargout > 1
        reg_grad = 2 * reg_weight * (pvec - opt.initialParameters);
    end
end

accum = opt.accumulateResiduals;
if isempty(accum)
    accum = struct('wells', [], 'types', [], 'steps', []);
end
nparam = cellfun(@(x)x.nParam, parameters);
p_org = pvec;
if opt.enforceBounds
    pvec = max(0, min(1, pvec));
end
pvec = mat2cell(pvec, nparam, 1);

% Create new setup, and set parameter values
pval = cell(size(parameters));
setupNew = setup;
setupNew.model.FlowDiscretization = [];
setupNew.model.FlowPropertyFunctions = [];
for k = 1:numel(parameters)
    pval{k}  = parameters{k}.unscale(pvec{k});
    setupNew = parameters{k}.setParameter(setupNew, pval{k});
end
[wellSols,states] = simulateScheduleAD(setupNew.state0, setupNew.model, setupNew.schedule,...
                                       'NonlinearSolver',opt.NonlinearSolver,...
                                       'Verbose',opt.Verbose, extra{:});
                                   
misfitVals = obj(setupNew.model, states, setupNew.schedule, states_ref, false, [],[]);
if ~isempty(accum.steps)
    tmp = repmat({0}, [max(accum.steps), 1]);
    for k = 1:numel(misfitVals)
        if accum.steps(k) > 0
            tmp{accum.steps(k)} = tmp{accum.steps(k)} + misfitVals{k};
        end
    end
    misfitVals = tmp;
end
misfitVals = (vertcat(misfitVals{:})).^(1/2);

% ADD REGULARIZATION TO THE OBJECTIVE VALUE
misfitVals_with_reg = misfitVals.^2 + reg_term;  % Since misfitVals are sqrt already

if nargout > 1
    objh = @(tstep,model,state) obj(setupNew.model, states, setupNew.schedule, states_ref, true, tstep, state);
    nms = applyFunction(@(x)x.name, parameters);
    scaledGradient = cell(numel(nms), 1);
    
    setupNew.model = setupNew.model.validateModel();
    gradient = computeSensitivitiesAdjointAD(setupNew, states, parameters, objh,...
                         'accumulateResiduals', accum, 'isScalar', false);
    for k = 1:numel(nms)
        scaledGradient{k} = parameters{k}.scaleGradient( gradient.(nms{k}), pval{k});
    end
    J = vertcat(scaledGradient{:})'/opt.objScaling;
    % adjust jacobian for sum r^2 -> sum r
    nzIx  = abs(misfitVals) > eps*norm(misfitVals);
    numnz = nnz(nzIx);
    J(nzIx,:) = spdiags((2*misfitVals(nzIx)).^(-1), 0, numnz, numnz)*J(nzIx,:);
    
    % ADD REGULARIZATION TO THE GRADIENT
    if ~isempty(reg_grad)
        J = J + reg_grad';
    end
    
    varargout{1} = J;
    if nargout > 2
        [varargout{2:3}] = deal(wellSols, states);
    end
    if nargout > 4
        varargout{4} =  setupNew;
    end
end

% RETURN THE REGULARIZED OBJECTIVE VALUE
misfitVals = misfitVals_with_reg;
end