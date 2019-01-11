function [G, data, Gs, valid_ix] = readAndPrepareForPostProcessorMRST(problem, steps, info, precomp)

model = problem.SimulatorSetup.model;
schedule = problem.SimulatorSetup.schedule;
G = model.G;
Gs = model.G;

init.PORO.values = model.rock.poro;
init.PERMX.values = model.rock.perm(:,1);
init.PERMY.values = model.rock.perm(:,2);
init.PERMZ.values = model.rock.perm(:,3);
init.DEPTH.values = G.cells.centroids(:,3);
init.PORV.values = G.cells.volumes;

% more fields from init might be interesting
data = setStatic([], init, {'PORO', 'PERMX', 'PERMY', 'PERMZ', 'NTG', 'DEPTH', 'PORV'});

% time in days (start and end of restart step)
startday = datenum(info.date(1, [3 2 1]));
data.time.prev = startday + info.time( max(steps-1,1) ) - info.time(1);
data.time.cur  = startday + info.time( steps ) - info.time(1);

[ws, states, report] = getPackedSimulatorOutput(problem);
t = cumsum(schedule.step.val);
for i = 1:numel(schedule.step.val)
    wells{i} = schedule.control(schedule.step.control(i)).W;
end




if isempty(precomp)
    data.states = states;
else
    data.states = cellfun(@(x)x.states{1}, precomp, 'UniformOutput', false);
end




valid_ix = isValidState(data.states);
if ~all(valid_ix)
    ns = nnz(~valid_ix);
    warning('Current version requires at least one open well.\n Skipping %2.0d of the %2.0d selected restart steps.', ns, numel(valid_ix));
    data.time.prev = data.time.prev(valid_ix);
    data.time.cur  = data.time.cur(valid_ix);
    data.states    = data.states(valid_ix);
end

% include some more fields from restart later on
data = setDynamic(G, data, valid_ix);

% set empty computed-field
data.computed = repmat(struct('name', '', 'values', [], 'limits', []), [0 1]);

end 

function data = setStatic(~, init, propnames)
for k = 1:numel(propnames)
    if isfield(init, propnames{k})
        v = init.(propnames{k}).values;
        data.static(k) = struct('name',    propnames{k}, ...
                                'values' , v, ...
                                'limits' , [min(v), max(v)]);
    end
end
end

function data = setDynamic(G, data, valid_ix)
if ~any(valid_ix)
    data.dynamic = [];
else
    flds = {'PRESSURE', 'SWAT', 'SOIL', 'SGAS'};
    p = cellfun(@(x)x.pressure, data.states, 'UniformOutput', false);
    p = horzcat(p{:})/barsa;
    data.dynamic(1) = struct('name', 'PRESSURE', ...
        'values' , p, ...
        'limits' , [min(min(p)), max(max(p))]);
    
    for k = 1: size(data.states{1}.s, 2)
        nm = flds{k+1};
        vals = cellfun(@(x)x.s(:,k), data.states, 'UniformOutput', false);
        vals = horzcat(vals{:});
        % take min/max over all steps
        [minv, maxv] = deal(min(min(vals)), max(max(vals)));
        data.dynamic(k+1) = struct('name', nm, ...
            'values' , vals, ...
            'limits' , [minv, maxv]);
    end
    
% Maybe include this in future for MRST input.
%     % additional fields of size corressponding to G.cells.num
%     ix = structfun(@(fld)all(cellfun(@numel, fld)==G.cells.num), rstrt);
%     fn = fieldnames(rstrt);
%     % dont include sats, pressure, flux
%     for pat = [flds, {'FLR', 'FLO'}]
%         ix = and(ix, ~strncmp(fn, pat, numel(pat{1})));
%     end
%     ix = find(ix);
%     k0 = numel(data.dynamic);
%     for k = 1:numel(ix)
%         nm = fn{ix(k)};
%         vals = horzcat(rstrt.(nm){:});
%         vals = vals(:, valid_ix);
%         [minv, maxv] = deal(min(min(vals)), max(max(vals)));
%         data.dynamic(k0+k) = struct('name', nm, ...
%             'values' , vals, ...
%             'limits' , [minv, maxv]);
%     end
end
end

function flag = isValidState(states)
flag = true(numel(states), 1);
for k = 1:numel(states)
    ws = states{k}.wellSol;
    stat = and([ws.status], abs([ws.val])>0);
    openPrd = any(and(stat, [ws.sign]<0));
    openInj = any(and(stat, [ws.sign]>0));
    flag(k) = openPrd || openInj;
end
end

function states = addConnectionPhaseFluxes(states, fluid, runspec)
% from output we typically only have total volume connection fluxes. 
% Approx water flux by using bW at bhp 

assert(runspec.OIL && runspec.WATER, 'Current code assumes both oil and water present')
for sk = 1:numel(states)
    ws = states{sk}.wellSol;
    for wk = 1:numel(ws)
        if size(ws(wk).flux, 2) == 1
            c = ws(wk).cells;
            resflux = ws(wk).flux;
            qwr = ws(wk).cqs(:,1)./states{sk}.b(c,1);
            qor = ws(wk).cqs(:,2)./states{sk}.b(c,2);
            if runspec.VAPOIL
                qor = qor - ws(wk).cqs(:,3).*states{sk}.rv(c)./states{sk}.b(c,2);
            end
            if ~runspec.GAS
                ws(wk).flux = [qwr, qor];
            else
                qgr = ws(wk).cqs(:,3)./states{sk}.b(c,3);
                if runspec.DISGAS
                    qgr = qgr - ws(wk).cqs(:,2).*states{sk}.rs(c)./states{sk}.b(c,3);
                    if runspec.DISGAS && runspec.VAPOIL
                        r = 1-states{sk}.rs(c).*states{sk}.rv(c);
                        [qor, qgr] = deal(qor./r, qgr./r);
                    end
                end
                ws(wk).flux = [qwr, qor, qgr];
            end
            % if remove crossflow ...
            %ws(wk).flux(resflux==0,:) = 0;
        end
    end
    states{sk}.wellSol = ws;
end
end

% function m = getModelDetails(init)
% [ih, lg] = deal(init.INTEHEAD.values, init.LOGIHEAD.values);
% unms = {'metric', 'field', 'lab'};
% m.units = unms{ih(3)};
% phnms = {'oil', 'water', 'gas'};
% phM   = [1 0 1 0 0 1 1
%          0 1 1 0 1 0 1
%          0 0 0 1 1 1 1];
% actPh = phM(:, ih(15));
% for k =1:3
%     m.(phnms{k}) = logical(actPh(k));
% end
% m.disgas = 
% m = 1;
% end
        
    

