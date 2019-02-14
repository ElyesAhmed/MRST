function [problem, state] = transportEquationBlackOilDG(state0, state, model, dt, drivingForces, varargin)

    opt = struct('Verbose'      , mrstVerbose, ...
                 'reverseMode'  , false      , ...
                 'scaling'      , []         , ...
                 'resOnly'      , false      , ...
                 'solveForWater', false      , ...
                 'solveForOil'  , true       , ...
                 'solveForGas'  , true       , ...
                 'iteration'    , -1         , ...
                 'stepOptions'  , []         ); % Compatibility only
    opt      = merge_options(opt, varargin{:});
    
    % Frequently used properties
    op       = model.operators;
    fluid    = model.fluid;
    rock     = model.rock;
    G        = model.G;
    disc     = model.disc;
    W        = drivingForces.W;
    disgas   = model.disgas;
    vapoil   = model.vapoil;
    psi      = disc.basis.psi;
    gradPsi  = disc.basis.grad_psi;
    
    % We may solve for both oil and water simultaneously
    solveAllPhases = opt.solveForWater && opt.solveForOil && opt.solveForGas;
    
    % Prepare state for simulation-----------------------------------------
    if opt.iteration == 1 && ~opt.resOnly 
        if model.tryMaxDegree
            % If we are at the first iteration, we try to solve using
            % maximum degree in all cells
            state.degree(~G.cells.ghost) = disc.degree;
        end
        if ~isempty(W)
            state.degree(vertcat(W.cells)) = 0;
        end
        % For cells that previously had less than nDof unknowns, we must
        % map old dofs to new
        state = disc.mapDofs(state, state0, 's');
        state = disc.mapDofs(state, state0, 'rs');
        state = disc.mapDofs(state, state0, 'rv');
        
    end
    % Update discretizaiton information. This is carried by the state
    % variable, and holds the number of dofs per cell + dof position in
    % state.sdof
    state0 = disc.updateDofPos(state0);
    [state, disc]  = disc.updateDofPos(state);
    %----------------------------------------------------------------------
    
    % Properties from current and previous timestep------------------------
    [p , sWdof , sOdof , sGdof , rSdof , rVdof , wellSol] = model.getProps(state , ...
                  'pressure', 'swdof', 'sodof', 'sgdof', 'rsdof', 'rvdof', 'wellsol');
    [p0, sWdof0, sOdof0, sGdof0, rSdof0, rVdof0         ] = model.getProps(state0, ...
                  'pressure', 'swdof', 'sodof', 'sgdof', 'rsdof', 'rvdof'           );
    % If timestep has been split relative to pressure, linearly interpolate
    % in pressure.
    if isfield(state, 'timestep')
        dt_frac = dt/state.timestep;
        p       = p.*dt_frac + p0.*(1-dt_frac);
    end
    %----------------------------------------------------------------------
    
    % Initialization of independent variables -----------------------------
    
    [sW , sO , sG ] = model.getProps(state , 'water', 'oil', 'gas');
    [sW0, sO0, sG0] = model.getProps(state0, 'water', 'oil', 'gas');
    st  = model.getCellStatusVO(state , sO ,  sW , sG );
    st0 = model.getCellStatusVO(state0, sO0,  sW0, sG0);
    
    if disgas || vapoil
        st  = cellfun(@(st) expandSingleValue(st, G), st , 'unif', false);
        st0 = cellfun(@(st) expandSingleValue(st, G), st0, 'unif', false);
        st  = cellfun(@(st) rldecode(st, state.nDof, 1), st, 'unif', false);
        st0 = cellfun(@(st) rldecode(st, state.nDof, 1), st0, 'unif', false);
        xDof = st{1}.*rSdof + st{2}.*rVdof + st{3}.*sGdof;
        gVar = 'xDof';
    else
        gVar = 'sGdof';
    end
    
    assert(~opt.reverseMode, 'Backwards solver not supported for splitting');
    if solveAllPhases
        if ~opt.resOnly
            if disgas || vapoil
                [sWdof, sOdof, xDof] = model.AutoDiffBackend.initVariablesAD(sWdof, sOdof, xDof);
            else
                [sWdof, sOdof, sGdof] = model.AutoDiffBackend.initVariablesAD(sWdof, sOdof, sGdof);
            end
        end
        primaryVars = {'sWdof', 'sOdof', gVar};
        sTdof = sOdof + sWdof + sGdof;
    else
        if ~opt.resOnly
            if disgas || vapoil
                [sWdof, xDof] = model.AutoDiffBackend.initVariablesAD(sWdof, xDof);
            else
                [sWdof, sGdof] = model.AutoDiffBackend.initVariablesAD(sWdof, sGdof);
            end
        end
        primaryVars = {'sWdof', gVar};
        sOdof     = -(sWdof + sGdof);
        ix        = disc.getDofIx(state, 1, Inf);
        sOdof(ix) = 1 + sOdof(ix);
        sTdof     = zeros(size(double(sWdof)));
        sTdof(ix) = 1;
    end
    
    if disgas || vapoil
        sGdof = st{2}.*(-sWdof);
        ix    = disc.getDofIx(state, 1, Inf);
        sGdof(ix) = sGdof(ix) + 1;
        sGdof = sGdof + st{3}.*xDof;
        if disgas
            rsSat = rldecode(fluid.rsSat(p), state.nDof, 1);
            rSdof = (~st{1}).*rsSat + st{1}.*xDof;
        end
        if vapoil
            rvSat = rldecode(fluid.rvSat(p), state.nDof, 1);
            rVdof = (~st{2}).*rvSat + st{2}.*xDof;
        end
    end
    disc.sample = sWdof;
    
    %----------------------------------------------------------------------

    % Pressure and saturation dependent properties-------------------------
    % Get multipliers
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    pvMult  = expandSingleValue(pvMult , G);
    pvMult0 = expandSingleValue(pvMult0, G);
    mobMult = expandSingleValue(mobMult, G);
    T       = op.T.*transMult;
    tm = ones(G.faces.num,1);
    tm(op.internalConn) = transMult;
    T_all   = model.operators.T_all.*tm;
    
    % Phase properties
    gdz = model.getGravityGradient();
    [b, mu, rho, mob] = getDerivedPropertyFunctionsBO(model, p, mobMult, st);
    bW   = b{1};   bO   = b{2};   bG   = b{3};
    rhoW = rho{1}; rhoO = rho{2}; rhoG = rho{3};
    mobW = mob{1}; mobO = mob{2}; mobG = mob{3};

    [sW, sO, sG] = disc.getCellMean(state, sWdof, sOdof, sGdof);
    
    % Gravity flux
    
    x_f = G.faces.centroids(disc.internalConn);
    c_l = disc.N(:,1);
    c_r = disc.N(:,2);
    [sW_l, sO_l, sG_l, rS_l, rV_l] = disc.evaluateDGVariable(x_f, c_l, state, sWdof, sOdof, sGdof, rSdof, rVdof);
    [sW_r, sO_r, sG_r, rS_r, rV_r] = disc.evaluateDGVariable(x_f, c_r, state, sWdof, sOdof, sGdof, rSdof, rVdof); 
    
    gW = (rhoW(c_l, sW_l) + rhoW(c_r, sW_r))/2.*gdz;
    if isfield(fluid, 'pcOW')
        gW = gW - op.Grad(fluid.pcOW(sW));
    end
    gO = (rhoO(c_l, rS_l) + rhoO(c_r, rS_r))/2.*gdz;
    gG = (rhoG(c_l, sG_l, rV_l) + rhoG(c_r, sG_r, rV_r))/2.*gdz;
    if isfield(fluid, 'pcOG')
        gG = gG + op.Grad(fluid.pcOG(sG));
    end
    P = sparse(find(op.internalConn), 1:nnz(op.internalConn), 1, G.faces.num, nnz(op.internalConn));
    gW = P*gW; gO = P*gO; gG = P*gG;
    % Add gravity flux where we have BCs to get correct cell values
    bc = drivingForces.bc;
    if ~isempty(bc)
        error('bcs not implemented');
        BCcells = sum(G.faces.neighbors(bc.face,:), 2);
        dz = G.cells.centroids(BCcells, :) - G.faces.centroids(bc.face,:);
        g = model.getGravityVector();
        rhoWBC = rhoW(BCcells);
        rhoOBC = rhoO(BCcells);
        rhoGBC = rhoO(BCcells);
        gW(bc.face) = rhoWBC.*(dz*g');
        gO(bc.face) = rhoOBC.*(dz*g');
    end    

    [qb_c, qb_f, r_c, r_fg] = computeSequentialFluxesDG(disc, model, state, T, T_all, ...
        {gW, gO, gG}, {mobW, mobO, mobG}, {bW, bO, bG}, {sWdof, sOdof, sGdof, sTdof}, {0*rSdof, rSdof, rVdof});
    [bWqW_c, bOqO_c, bGqG_c] = deal(qb_c{:});
    [bWqW_f, bOqO_f, bGqG_f] = deal(qb_f{:});
    % Well contributions---------------------------------------------------
    if ~isempty(W)
        % Total well flux, composition and mappings
        perf2well = getPerforationToWellMapping(W);
        wc        = vertcat(W.cells);
        wflux     = zeros(G.cells.num,1);
        wflux(wc) = sum(vertcat(wellSol.flux), 2)./G.cells.volumes(wc);
        isInj     = wflux > 0;
        compWell  = vertcat(W.compi);
        compPerf  = zeros(G.cells.num, 3);
        compPerf(wc,:) = compWell(perf2well,:);
        
        % Saturations at cubature points
        [~, x_w, c_w] = disc.getCubature(wc, 'volume');
        [sW_w, sO_w, sG_w, sT_w, rS_w, rV_w] = disc.evaluateDGVariable(x_w, c_w, state, sWdof, sOdof, sGdof, sTdof, rSdof, rVdof);
        
        mobW_w = mobW(c_w, sW_w, sT_w);
        mobO_w = mobO(c_w, sO_w, sT_w, rS_w);
        mobG_w = mobG(c_w, sG_w, sT_w, rV_w);
        mobT_w = mobW_w + mobO_w + mobG_w;
        
        fW_w = ~isInj(c_w).*sT_w.*mobW_w./mobT_w + isInj(c_w).*compPerf(c_w,1);
        fO_w = ~isInj(c_w).*sT_w.*mobO_w./mobT_w + isInj(c_w).*compPerf(c_w,2);
        fG_w = ~isInj(c_w).*sT_w.*mobG_w./mobT_w + isInj(c_w).*compPerf(c_w,3);
        
        bWqW_w = bW(c_w, sW_w).*wflux(c_w).*sT_w.*fW_w;
        bOqO_w = bO(c_w, rS_w).*wflux(c_w).*sT_w.*fO_w;
        bGqG_w = bG(c_w, sG_w, rV_w).*wflux(c_w).*sT_w.*fG_w;
        
        % Water well contributions
        integrandW = @(psi, gradPsi) bWqW_w.*psi;
        srcW_w = disc.cellInt(integrandW, wc, state, sWdof);
        
        % Oil well contributions
        integrandO = @(psi, gradPsi) (bOqO_w + rV_w.*bGqG_w).*psi;
        srcO_w = disc.cellInt(integrandO, wc, state, sWdof);
        
        % Oil well contributions
        integrandG = @(psi, gradPsi) (bGqG_w + rS_w.*bOqO_w).*psi;
        srcG_w = disc.cellInt(integrandG, wc, state, sWdof);
                                 
        % Store well fluxes
        ix     = disc.getDofIx(state, 1, wc);
        wfluxW = double(srcW_w(ix));
        wfluxO = double(srcO_w(ix));
        wfluxG = double(srcG_w(ix));
        for wNo = 1:numel(W)
            perfind = perf2well == wNo;
            state.wellSol(wNo).qWs = sum(wfluxW(perfind));
            state.wellSol(wNo).qOs = sum(wfluxO(perfind));
            state.wellSol(wNo).qGs = sum(wfluxG(perfind));
        end

    end
    %----------------------------------------------------------------------

    % Evaluate saturation at cubature points-------------------------------
    % Cell cubature points
    [~, x_c, c] = disc.getCubature((1:G.cells.num)', 'volume');
    [sW_c , sO_c , sG_c , rV_c , rS_c] = disc.evaluateDGVariable(x_c, c, state , sWdof , sOdof , sGdof, rVdof , rSdof);
    [sW0_c, sO0_c, sG0_c, rV0_c, rS0_c]       = disc.evaluateDGVariable(x_c, c, state0, sWdof0, sOdof0, sGdof0, rVdof0, rSdof0);
    % Face cubature points
    [~, x_f, ~, f] = disc.getCubature((1:G.cells.num)', 'surface'); 
    %----------------------------------------------------------------------
    
    [eqs, names, types] = deal(cell(1,2 + solveAllPhases));
    [types{:}] = deal('cell');
    
    eqNo = 1;
    % Water equation-------------------------------------------------------
    if opt.solveForWater

        mW  = pvMult(c).*rock.poro(c).*bW(c, sW_c).*sW_c;
        mW0 = pvMult(c).*rock.poro(c).*bW(c, sW0_c).*sW0_c;

        water =   disc.inner((mW - mW0)/dt, psi, 'dV') ...
                - disc.inner(bWqW_c, gradPsi, 'dV')    ...
                + disc.inner(bWqW_f, psi, 'dS');

        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            water(ix) = water(ix) - srcW_w(ix);
        end
        eqs{eqNo}   = water;
        names{eqNo} = 'water';
        eqNo        = eqNo + 1;
    end
    %----------------------------------------------------------------------
    
    % Oil equation---------------------------------------------------------
    if opt.solveForOil
        % Cell values
                
        mO  = pvMult(c).*rock.poro(c).*bO(c, sO_c, rS_c).*sO_c;
        mO0 = pvMult(c).*rock.poro(c).*bO(c, sO0_c, rS0_c).*sO0_c;
        
        oil  =   disc.inner((mO - mO0)/dt, psi, 'dV') ...
              - disc.inner(bOqO_c, gradPsi, 'dV')    ...
              + disc.inner(bOqO_f, psi, 'dS');
        
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            oil(ix) = oil(ix) - srcO_w(ix);
        end
        eqs{eqNo}   = oil;
        names{eqNo} = 'oil';
        eqNo        = eqNo + 1;
    end
    %----------------------------------------------------------------------

    % Gas equation---------------------------------------------------------
    if opt.solveForGas
        % Cell values
        
        mG  = pvMult(c).*rock.poro(c).*bG(c, sG_c, rV_c).*sG_c;
        mG0 = pvMult(c).*rock.poro(c).*bG(c, sG0_c, rV0_c).*sG0_c;
        
        gas =   disc.inner((mG - mG0)/dt, psi, 'dV') ...
              - disc.inner(bGqG_c, gradPsi, 'dV')    ...
              + disc.inner(bGqG_f, psi, 'dS');
        
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            gas(ix) = gas(ix) - srcG_w(ix);
        end
        eqs{eqNo}   = gas;
        names{eqNo} = 'gas';
    end
    %----------------------------------------------------------------------
    
    % Add BCs--------------------------------------------------------------
    if ~isempty(bc)
        % Boundary faces
        faces = bc.face;
        % Saturation outside boundary
        sL  = bc.sat;
        % Face cubature foordinates
        [~, x, ~, fBC] = disc.getCubature(faces, 'face');
        cBC = sum(G.faces.neighbors(fBC,:),2);
%         [xR, ~, ~] = disc.transformCoords(x, cBC);
        % Mapping from BC faces to global faces
        globFace2BCface        = nan(G.faces.num,1);
        globFace2BCface(faces) = 1:numel(faces);        
        locFaceNo = globFace2BCface(fBC);
        % Determine injectng boundaries
        sgn = 1 - 2*(G.faces.neighbors(faces, 1) == 0);
        isInj = vT(faces) > 0 & sgn < 0;
        % Upstream saturation
        [sW_r, sO_r, sTR]  = disc.evaluateDGVariable(x, cBC, state, sWdof, sOdof, sTdof);
        sWBC = sL(locFaceNo,1).*isInj(locFaceNo) + sW_r.*(~isInj(locFaceNo));
        sOBC = sL(locFaceNo,2).*isInj(locFaceNo) + sO_r.*(~isInj(locFaceNo));
        sTBC = sum(sL(locFaceNo,:),2).*isInj(locFaceNo) + sTR.*(~isInj(locFaceNo));
        % Frational flow functions
        fWBC = fW(sWBC, sOBC, sTBC, cBC, cBC);
        fOBC = fO(sWBC, sOBC, sTBC, cBC, cBC);
        if opt.solveForWater
            % Add water bc flux to water equation
            faceIntegrand = @(psi) (bW(cBC).*fWBC.*vT(fBC) ...
                      + bW(cBC).*fWBC.*mobO(sOBC,sTBC,cBC).*(TgW(fBC) - TgO(fBC))).*psi;
            fluxWBC = disc.faceFluxIntBC(faceIntegrand, bc, state, sWdof);
            water   = water + fluxWBC;
        end
        if opt.solveForOil
            % Add oil bc flux to oil equation
            faceIntegrand = @(psi) (bO(cBC).*fOBC.*vT(fBC) ...
                      + bO(cBC).*fOBC.*mobW(sWBC,sTBC,cBC).*(TgO(fBC) - TgW(fBC))).*psi;
            fluxOBC = disc.faceFluxIntBC(faceIntegrand, bc, state, sOdof);
            oil     = oil + fluxOBC;
        end
    end
    %----------------------------------------------------------------------
    
    % Add sources----------------------------------------------------------
    src = drivingForces.src;
    if ~isempty(src) 
        % Cubature
        [~, ~, cSRC] = disc.getCubature(src.cell, 'volume');
        % Mapping from source cells to global cells
        globCell2SRCcell = nan(G.cells.num,1);
        globCell2SRCcell(src.cell) = 1:numel(src.cell);
        cSRCloc = globCell2SRCcell(cSRC);
        % Total rate and saturaion at cubature points
        qT   = src.rate(cSRCloc)./G.cells.volumes(cSRC);
        sSRC = src.sat(cSRCloc,:);
        if opt.solveForWater
            % Add water source to water equation
            srcIntegrand = @(psi, gradPsi) bW(cSRC).*qT.*sSRC(:,1).*psi;
            srcW  = disc.cellInt(srcIntegrand, src.cell, state, sWdof);
            water = water - srcW;
        end
        if opt.solveForOil
            % Add oil source to oil equation
            srcIntegrand = @(psi, gradPsi) bO(cSRC).*qT.*sSRC(:,2).*psi;
            srcO = disc.cellInt(srcIntegrand, src.cell, state, sOdof);
            oil  = oil - srcO;
        end
    end
    %----------------------------------------------------------------------
    
    % Make Linearized problem----------------------------------------------
    % Scale equations
    if ~model.useCNVConvergence
        pv = rldecode(op.pv, state.nDof, 1);
        for eqNo = 1:numel(eqs)
            eqs{eqNo} = eqs{eqNo}.*(dt./pv);
        end
    end
    % Extra state output
    if model.extraStateOutput
        state     = model.storeDensity(state, rhoW, rhoO, []);
        state.cfl = dt.*sum(abs(vTc)./G.cells.dx,2);
    end
    % Linearize
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
    %----------------------------------------------------------------------
    
    % Extract subproblem if we are solving subproblem----------------------
    if any(strcmpi(G.type, 'subgrid'))
        ix = disc.getDofIx(state, Inf, ~G.cells.ghost);
        
        for eqNo = 1:numel(problem.equations)
            eq = problem.equations{eqNo};
            if isa(eq, 'ADI')
                eq.val = eq.val(ix);
                eq.jac = cellfun(@(j) j(ix,ix), eq.jac, 'unif', false);
            else
                eq = eq(ix);
            end
            problem.equations{eqNo} = eq;
        end
    end
    %----------------------------------------------------------------------
    
end

% Expang single scalar values to one per cell------------------------------
function v = expandSingleValue(v,G)
    if numel(double(v)) == 1
        v = v*ones(G.cells.num,1);
    end
end
%--------------------------------------------------------------------------