function [problem, state] = equationsOilWater(state0, state, model, dt, drivingForces, varargin)
% Get linearized problem for oil/water system with black oil-style
% properties

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', [],...
             'iteration', -1, ...
             'stepOptions', [],...
             'addflux',false);  % Compatibility only

opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

W = drivingForces.Wells;
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

s = model.operators;
G = model.G;
f = model.fluid;

hst = opt.history;

[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');

[p0, sW0] = model.getProps(state0, 'pressure', 'water');


pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, qWs, qOs, pBH] = ...
            initVariablesADI(p, sW, qWs, qOs, pBH);
    else
        [p0, sW0, tmp, tmp, tmp] = ...
            initVariablesADI(p0, sW0,          ...
            zeros(size(qWs)), ...
            zeros(size(qOs)), ...
            zeros(size(pBH)));                          %#ok
    end
end
primaryVars = {'pressure', 'sW', 'qWs', 'qOs', 'bhp'};

clear tmp
g  = norm(gravity);

%--------------------
%check for p-dependent tran mult:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end
transMult=1;
if isfield(f, 'transMult')
   transMult=f.transMult(p); 
end
%check for capillary pressure (p_cow)
pcOW = 0;
if isfield(f, 'pcOW')
    pcOW  = f.pcOW(sW);
end

trans=s.T.*transMult;
% -------------------------------------------------------------------------
[krW, krO] = f.relPerm(sW);
%krW = f.krW(sW);
%krO = f.krO(1-sW);

% water props (calculated at oil pressure OK?)
%bW     = f.bW(p);
bW     = f.bW(p-pcOW);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
%mobW   = trMult.*krW./f.muW(p);
mobW   = trMult.*krW./f.muW(p-pcOW);
dpW     = s.grad(p-pcOW) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
% water upstream-index
upc = (double(dpW)>=0);
bWvW = s.faceUpstr(upc, bW.*mobW).*trans.*dpW;


% oil props
bO     = f.bO(p);
rhoO   = bO.*f.rhoOS;
rhoOf  = s.faceAvg(rhoO);
dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
% oil upstream-index
upc = (double(dpO)>=0);
if isfield(f, 'BOxmuO')
    % mob0 is already multplied with b0
    mobO   = trMult.*krO./f.BOxmuO(p);
    bOvO   = s.faceUpstr(upc, mobO).*trans.*dpO;
else
    mobO   = trMult.*krO./f.muO(p);
    bOvO   = s.faceUpstr(upc, bO.*mobO).*trans.*dpO;
end
if(opt.addflux)
   state.bOvO=bOvO;
   if isfield(f, 'BOxmuO')
     state.OvO=s.faceUpstr(upc, mobO./b0).*trans.*dpO;
   else
     state.OvO=s.faceUpstr(upc, mobO).*trans.*dpO;
   end
   state.bWvW=bWvW;
   state.WvW= s.faceUpstr(upc, mobW).*trans.*dpW;
end

% EQUATIONS ---------------------------------------------------------------
% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sW) - pvMult0.*f.bO(p0).*(1-sW0) ) - s.div(bOvO);

% water:
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0 ) - s.div(bWvW);

names = {'oil', 'water'};
types = {'cell', 'cell'};
% well equations
if ~isempty(W)
    if ~opt.reverseMode
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS, f.rhoOS];
        bw   = {bW(wc), bO(wc)};
        mw   = {mobW(wc), mobO(wc)};
        s = {sW, 1 - sW};

        wm = WellModel();
        [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, W, wellSol, ...
                                             pBH, {qWs, qOs}, pw, rhos, bw, mw, s, {},...
                                             'nonlinearIteration', opt.iteration);
        eqs(3:4) = weqs;
        eqs{5} = ctrleqs;
        
        eqs{1}(wc) = eqs{1}(wc) - cqs{2};
        eqs{2}(wc) = eqs{2}(wc) - cqs{1};
        
        names(3:5) = {'oilWells', 'waterWells', 'closureWells'};
        types(3:5) = {'perf', 'perf', 'well'};
    else
        % in reverse mode just gather zero-eqs of correct size
        for eqn = 3:5
            nw = numel(state0.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(3:5) = {zw, zw, zw};
        end
        names(3:5) = {'empty', 'empty', 'empty'};
        types(3:5) = {'none', 'none', 'none'};
    end
else % no wells
    eqs(3:5) = {pBH, pBH, pBH};  % empty  ADIs
    names(3:5) = {'empty', 'empty', 'empty'};
    types(3:5) = {'none', 'none', 'none'};
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
problem.iterationNo = opt.iteration;
end
%--------------------------------------------------------------------------









