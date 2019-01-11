function [] = processStatesDiagnostics(problem, varargin)
opt = struct('outputdir',         [], ...
             'multiple',        true, ...
             'maxTOF',     500*year, ...
             'startdate',    [0 0 0]);

datadir = problem.OutputHandlers.states.getDataPath;

if isempty(opt.outputdir)
    opt.outputdir = fullfile(datadir, 'mrst_diagnostics');
end

if exist(opt.outputdir,'dir')~=7
    flag = mkdir(opt.outputdir);
    if ~flag
        error('Unable to create directory...')
    end
end

prefx = problem.BaseName;


%[fluid, pvtdeck] = initFluidFromOutput(init);

% only simulation-grid here
G = problem.SimulatorSetup.model.G;
schedule = problem.SimulatorSetup.schedule;


[ws, states, report] = getPackedSimulatorOutput(problem);
t = cumsum(schedule.step.val)./day;
for i = 1:numel(schedule.step.val)
    wells{i} = schedule.control(schedule.step.control(i)).W;
end

startdate = opt.startdate;
steps = 1:numel(t);
startday = datenum(startdate(1, [3 2 1]));
time.prev = startday + t( max(steps-1,1) ) - t(1);
time.cur  = startday + t( steps ) - t(1);

% startday = datenum(startdate);
% time.prev = [startday; startday; startday + (schedule.step.val(1:end-1)./day)];
% time.cur  = startday + schedule.step.val(:)./day;


if opt.multiple
    t0 = tic;
    fprintf('Computing diagnostics:     ')
    fnm = @(k)fullfile(opt.outputdir, sprintf([prefx,'_diagn%0.4d.mat'], k));
    for k = 1:numel(states)
        fprintf('\b\b\b\b%3.1d%%', round(100*k/numel(states)));
        tmptime = struct('cur', time.cur(k), 'prev', time.prev(k));
        cur = struct('states', {states(k)}, 'wells', {wells(k)}, 'diagnostics', [], 'time', tmptime);
        % switch off verbose here
        vb = mrstVerbose;
        mrstVerbose('off');
        cur = computeDiagnostics(G, cur, opt.maxTOF); %#ok
        mrstVerbose(vb);
        save(fnm(k), '-struct', 'cur');
        dt = toc(t0);
        time_left =  (numel(states)-k+1)*dt/k; %#ok
    end
    fprintf(', done\n')
else
    error('currently supports only multiple')
end
end


function states = addConnectionPhaseFluxes(states, fluid, runspec)
% from output we typically only have total volume connection fluxes. 
% Approx water flux by using bW at bhp 
if runspec.DISGAS && runspec.VAPOIL
    warning('Both occuring DISGAS and VAPOIL, connection fluxes will not be accurate')
end
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
    

