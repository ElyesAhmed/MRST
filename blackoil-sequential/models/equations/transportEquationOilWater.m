function [problem, state] = transportEquationOilWater(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', [],...
             'iteration', -1, ...
             'stepOptions', []);  % Compatibility only

opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

s = model.operators;
f = model.fluid;
G = model.G;

[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');

[p0, sW0] = model.getProps(state0, 'pressure', 'water');

wflux = vertcat(wellSol.flux);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        sW = initVariablesADI(sW);
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'sW'};

clear tmp
g  = norm(gravity);


% -------------------------------------------------------------------------
[krW, krO] = f.relPerm(sW);

clear krW_o krO_o

dZ = s.grad(G.cells.centroids(:,3));
g = norm(gravity);

sO = 1 - sW;
% Water
[bW, rhoW, mobW, dpW, Gw] = propsOW_water(sW, krW, g, dZ, f, p, s);
[bO, rhoO, mobO, dpO, Go] = propsOW_oil(  sO, krO, g, dZ, f, p, s);

if 1
    upco = state.upstream(:, 2);
    upcw = state.upstream(:, 1);
else
    upcw = (double(dpW)>=0);
    upco = (double(dpO)>=0);
end

intx = ~any(G.faces.neighbors == 0, 2);
vT = state.flux(intx);

mobOf = s.faceUpstr(upco, mobO);
mobWf = s.faceUpstr(upcw, mobW);

f_w = mobWf./(mobOf + mobWf);
bWvW   = s.faceUpstr(upcw, bW).*f_w.*(vT + s.T.*mobOf.*(Go - Gw));


if ~isempty(W)
    perf2well = getPerforationToWellMapping(W);
    wc = vertcat(W.cells);
    
    mobWw = mobW(wc);
    mobOw = mobO(wc);
    totMobw = mobWw + mobOw;

    f_w_w = mobWw./totMobw;
    f_o_w = mobOw./totMobw;

    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    f_w_w(isInj) = compPerf(isInj, 1);
    f_o_w(isInj) = compPerf(isInj, 2);

    bWqW = bW(wc).*f_w_w.*wflux;
    bOqO = bO(wc).*f_o_w.*wflux;

    % Store well fluxes
    wflux_O = double(bOqO);
    wflux_W = double(bWqW);
    
    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(wflux_O(perfind));
        state.wellSol(i).qWs = sum(wflux_W(perfind));
    end

end
%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end
wat = (s.pv/dt).*(pvMult.*bW.*sW       - pvMult0.*f.bW(p0).*sW0    ) - s.div(bWvW);
wat(wc) = wat(wc) - bWqW;

eqs{1} = wat;
types = {'cell'};
names = {'water'};

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
problem.iterationNo = opt.iteration;

% perf2well = getPerforationToWellMapping(W);

end
