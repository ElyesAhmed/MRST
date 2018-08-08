classdef DGDiscretization < HyperbolicDiscretization
    
    properties

        dim
        degree
        basis
        volumeCubature
        surfaceCubature
        useUnstructCubature
        jumpTolerance
        outTolerance
        meanTolerance
        internalConnParent
        limiterType
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function disc = DGDiscretization(model, dim, varargin)
            
            disc = disc@HyperbolicDiscretization(model);
            
            disc.dim    = dim;
            disc.degree = 1;
            disc.basis  = 'legendre';
            disc.jumpTolerance = 0.2;
            disc.outTolerance = 1e-4;
            disc.meanTolerance = 1e-4;
            disc.useUnstructCubature = false;
            disc.internalConnParent = disc.internalConn;
            disc.limiterType = 'kill';
            
            disc        = merge_options(disc, varargin{:});
            
            disc.basis  = dgBasis(dim, disc.degree, disc.basis);
            
            G = disc.G;
            disc.G.parent = G;
            if G.griddim == 2
                if disc.degree == 0 || disc.useUnstructCubature
                    disc.volumeCubature = Unstruct2DCubature(disc.G, disc.degree + 1, disc.internalConn);
                else
                    disc.volumeCubature  = TriangleCubature(disc.G, disc.degree + 1, disc.internalConn);
                end
                disc.surfaceCubature  = LineCubature(disc.G, disc.degree + 1, disc.internalConn);
            else
                if disc.degree == 0 || disc.useUnstructCubature
                    disc.volumeCubature = Unstruct3DCubature(disc.G, disc.degree + 1, disc.internalConn);
                    disc.surfaceCubature = Unstruct2DCubature(disc.G, disc.degree + 1, disc.internalConn);
%                     disc.surfaceCubature = TriangleCubature(disc.G, disc.degree + 1, disc.internalConn);
                else
                    disc.volumeCubature  = TetrahedronCubature(disc.G, disc.degree + 1, disc.internalConn);
                    disc.surfaceCubature = TriangleCubature(disc.G, disc.degree + 1, disc.internalConn);
                end
            end
%             [W, x, w, ii, jj, cellNo] = disc.volumeCubature.getCubature((1:G.cells.num)', 'cell');
%             disc.volumeCubature.W = W;
%             [W, x, w, ii, jj, cellNo] = disc.surfaceCubature.getCubature((1:G.cells.num)', 'cellsurface');
%             disc.surfaceCubature.W = W;
            
        end
        
        %-----------------------------------------------------------------%        
        function state = assignDofFromState(disc, state)

            state.degree = repmat(disc.degree, disc.G.cells.num, 1);
            state = disc.updateDofPos(state);
            
            state.nDof = disc.getnDof(state);
            sdof = zeros(sum(state.nDof), size(state.s,2));
            
            ix = disc.getDofIx(state, 1);
            sdof(ix, :) = state.s;

            state.sdof = sdof;

        end
        
        %-----------------------------------------------------------------%
        function state = updateDofPos(disc, state)

            dp = reshape((1:disc.G.cells.num*disc.basis.nDof)', disc.basis.nDof, []);            
            if nargin == 1
                nd = repmat(disc.basis.nDof, numel(state.cells), 1);
            else
                nd = disc.getnDof(state);
                subt = cumsum([0; disc.basis.nDof - nd(1:end-1)]);
                [ii, jj, v] = find(dp);
                
                if size(ii,1) == 1, ii = ii'; end
                if size(jj,1) == 1, jj = jj'; end
                if size(v ,1) == 1, v  =  v'; end

                cnDof = cumsum(nd);

                v = v - subt(jj);
                v(v > cnDof(jj)) = 0;
                dp = full(sparse(ii, jj, v));
            end
            
            state.nDof   = nd;
            state.dofPos = dp;
            
%             disc.state = state;
            
        end
        
        %-----------------------------------------------------------------%
        function [xhat, translation, scaling] = transformCoords(disc, x, cells, inverse, useParent)
            
            G = disc.G;
            
            if nargin < 4, inverse = false; end
            if nargin < 5, useParent  = false; end
            
            if isfield(G, 'mappings') && useParent
%                 cells = G.mappings.cellMap.new2old(cells);
                G = G.parent;
            end
            
            translation = -G.cells.centroids(cells,:);
            if isfield(G.cells, 'dx')
                scaling = 1./(G.cells.dx(cells,:)/2);
            else
                scaling = 1./(G.cells.diameters(cells)/(2*sqrt(G.griddim)));
            end
            
            if ~inverse%nargin < 4 || ~inverse
                xhat = (x + translation).*scaling;
                xhat = xhat(:, 1:disc.dim);
                scaling = scaling(:, 1:disc.dim);
                translation = translation(:, 1:disc.dim);
            else
                xhat = x./scaling - translation;
            end
               
        end
        
        %-----------------------------------------------------------------%
        function v = trimValues(disc, v)
            
            tol = eps(mean(disc.G.cells.volumes));
            tol = -inf;
%             tol = 1e-7;
            ix = abs(v) < tol;
            if isa(v, 'ADI')
                v.val(ix) = 0;
            else
                v(ix) = 0;
            end
            
        end
        
        %-----------------------------------------------------------------%
        function ix = getDofIx(disc, state, dofNo, cells, includezero)

            G = disc.G;
            if nargin == 2
                cells = 1:disc.G.cells.num;
                dofNo = 1:disc.basis.nDof;
            elseif nargin == 3
                cells = 1:G.cells.num;
            else
                if isempty(dofNo)
                    dofNo = 1:disc.basis.nDof;
                end
            end
            
            ix = state.dofPos(dofNo, cells);
            ix = ix(:);
            
            if nargin < 5
                includezero = false;
            end
            if ~includezero
                ix(ix == 0) = [];
            end
              
        end
        
        %-----------------------------------------------------------------%
        function nDof = getnDof(disc, state)
            
            if disc.degree < 0
                nDof = 0;
            else
                nDof = factorial(state.degree + disc.dim)...
                       ./(factorial(disc.dim).*factorial(state.degree));
            end
            
        end
        
        %-----------------------------------------------------------------%
        function state = mapDofs(disc, state, state0)
            
            state = disc.updateDofPos(state);
            
            if all(state.nDof == state0.nDof)
                return
            else
                sdof = zeros(sum(state.nDof), size(state0.sdof,2));
                for dofNo = 1:disc.basis.nDof
                    ix      = disc.getDofIx(state, dofNo, (1:disc.G.cells.num)', true);
                    ix0 = disc.getDofIx(state0, dofNo, (1:disc.G.cells.num)', true);
                    
%                     sdof(ix(ix0 > 0),:) = state0.sdof(ix0(ix0 > 0),:);
                    sdof(ix(ix0 > 0 & ix > 0),:) = state0.sdof(ix0(ix0 > 0 & ix > 0),:);
                    
                end
                state.sdof = sdof;
            end

        end
        
        %-----------------------------------------------------------------%
        function s = evaluateSaturation(disc, x, cells, dof, state)
            
            psi     = disc.basis.psi;
            nDof    = state.nDof;
            nDofMax = disc.basis.nDof;
            
            ix = disc.getDofIx(state, 1, cells);
            s = dof(ix).*0;
            for dofNo = 1:nDofMax
                keep = nDof(cells) >= dofNo;
                ix = disc.getDofIx(state, dofNo, cells(keep));
                if all(keep)
                    s = s + dof(ix).*psi{dofNo}(x(keep,:));
                else
                    s(keep) = s(keep) + dof(ix).*psi{dofNo}(x(keep,:));
                end

            end
            
        end
        
        %-----------------------------------------------------------------%
        function state = getCellSaturation(disc, state)
            
            
%             [W, x, cellNo, faceNo] = disc.getCubature((1:disc.G.cells.num)', 'volume');
            [W, x, cellNo, faceNo] = disc.getCubature((1:disc.G.cells.num)', 'volume');
%             getCubature(disc, cells, type)
%             [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(disc.G, (1:disc.G.cells.num)', max(disc.degree), 'volume');
%             W = sparse(ii, jj, w);

            x = disc.transformCoords(x, cellNo);
            
            sdof = state.sdof;
            nPh = size(sdof,2);
            s = zeros(disc.G.cells.num, nPh);
            for phNo = 1:nPh
                s(:,phNo) = (W*disc.evaluateSaturation(x, cellNo, sdof(:,phNo), state))./disc.G.cells.volumes;
            end
            
            state.s = s;
            
        end
        
        %-----------------------------------------------------------------%
        function I = cellInt(disc, integrand, f, cells, sdof, sdof0, state, state0)
        
            G    = disc.G;
            psi  = disc.basis.psi;
            grad_psi = disc.basis.grad_psi;
            nDof = state.nDof;
            nDofMax = disc.basis.nDof;
            
            [W, x, cellNo, faceNo] = disc.getCubature(cells, 'volume');
            [x, ~, scaling] = disc.transformCoords(x, cellNo);

            s  = disc.evaluateSaturation(x, cellNo , sdof , state );
            s0 = disc.evaluateSaturation(x, cellNo, sdof0, state0);
            f = f(s, 1-s, cellNo, cellNo);
            
            I = integrand(sdof, sdof, sdof, ones(numel(double(sdof)), 1), 1, ones(1, disc.dim)).*0;
            
            for dofNo = 1:nDofMax
                
                keepCells = nDof(cells) >= dofNo;
                
                if any(keepCells)
                
                    ix = disc.getDofIx(state, dofNo, cells(keepCells));
                    i  = W*integrand(s, s0, f, cellNo, psi{dofNo}(x), grad_psi{dofNo}(x).*scaling);
                    I(ix) = i(keepCells);
                    
                elseif numel(cells) == disc.G.cells.num
                    
                    warning('No cells with %d dofs', dofNo);
                    
                end
                
            end
            
            I = disc.trimValues(I);
            
        end
        
        %-----------------------------------------------------------------%
        function I = faceFluxInt(disc, integrand, f, cells, sdof, state, T, vT, g, mob)
            
            psi     = disc.basis.psi;
            nDof    = state.nDof;
            nDofMax = disc.basis.nDof;
            
            [W, x, cellNo, faceNo] = disc.getCubature(cells, 'surface');
            nPh = numel(g);
            
            if isempty(faceNo)
                I = 0;
                return
            end
            
            [flag_v, flag_G, upCells_v, upCells_G, s_v, s_G] = disc.getUpstreamCell(faceNo, x, T, vT, g, mob, sdof, state);
            
%             [upCells_v, upCells_G] = deal(repmat({cR}, 1, nPh));
%             [x_v, s_v, x_G, s_G] = deal(cell(1, nPh));
%             for phNo = 1:nPh
%                 
%                 upCells_v{phNo}(flag_v(:,phNo)) = cL(flag_v(:,phNo));
%                 x_v{phNo} = disc.transformCoords(x, upCells_v{phNo});
%                 s_v{phNo} = disc.evaluateSaturation(x_v{phNo}, upCells_v{phNo}, sdof, state);
%                 
%                 upCells_G{phNo}(flag_G(:,phNo)) = cL(flag_G(:,phNo));
%                 x_G{phNo} = disc.transformCoords(x, upCells_G{phNo});
%                 s_G{phNo} = disc.evaluateSaturation(x_G{phNo}, upCells_G{phNo}, sdof, state);
%                 
%             end
            
            f_v = f(s_v{1}, 1-s_v{2}, upCells_v{1}, upCells_v{2});
            f_G = f(s_G{1}, 1-s_G{2}, upCells_G{1}, upCells_G{2});
            
            [x_c, ~, ~] = disc.transformCoords(x, cellNo);

            I = integrand(sdof, sdof, sdof, sdof, 1, 1, 1, 1, 1).*0;
            
            for dofNo = 1:nDofMax
                
                keepCells = nDof(cells) >= dofNo;
                
                if any(keepCells)
                    
                    ix = disc.getDofIx(state, dofNo, cells(keepCells)');
                    i  = W*integrand(f_v, f_G, s_v{1}, s_G{1}, cellNo, upCells_v{1}, upCells_G{1}, faceNo, psi{dofNo}(x_c));
                    I(ix) = i(keepCells);
                    
                elseif numel(cells) == disc.G.cells.num
                    
                    warning('No cells with %d dofs', dofNo);
                    
                end
                
            end
            
            I = disc.trimValues(I);
            
        end
        
        %-----------------------------------------------------------------%
        function I = faceFluxIntBC(disc, integrand, f, sdof, state, bc, mob, g, vT, T)
            
            G       = disc.G;
            psi     = disc.basis.psi;
            nDof    = state.nDof;
            nDofMax = disc.basis.nDof;
            
            isFlux = strcmpi(bc.type, 'flux');
            faces = bc.face(isFlux);
            
            faceMap = nan(G.faces.num,1);
            faceMap(bc.face) = 1:numel(faces);
            
            cellMap = nan(G.parent.cells.num,1);
            cellMap(bc.cell(:,1)) = 1:numel(bc.face);
            
            [W, x, cellNo, faceNo] = disc.getCubature(faces, 'face');
            
            cL = bc.cell(faceMap(faceNo),1);
            cR = bc.cell(faceMap(faceNo),2);
            xL = disc.transformCoords(x, cL, false, true);
            xR = disc.transformCoords(x, cR, false, true);
            
            cLloc = cL;
            cRloc = cR;
            if isfield(G, 'mappings')
                cLloc = cellMap(cL);
                cRloc = G.mappings.cellMap.old2new(cR);
            end
            
            sL = disc.evaluateSaturation(xL, cLloc, bc.sat(:,1), bc.state);
            sR = disc.evaluateSaturation(xR, cRloc, sdof, state);
            
            
            mm1 = [bc.mob(cellMap(cL),1); double(mob{1}(sR  , cRloc))];
            mm2 = [bc.mob(cellMap(cL),2); double(mob{2}(1-sR, cRloc))];
            
            mob{1} = mm1;
            mob{2} = mm2;
            
            g = cellfun(@(g) g(faceNo), g, 'unif', false);
            
            vT = -bc.value(faceMap(faceNo));
            T = T(faceNo);
            
            N = [1:numel(cL); numel(cL)+1:2*numel(cL)]';
            upw = @(flag, x)faceUpstr(flag, x, N, [size(N,1), max(max(N))]);
            [flag_v, flag_G] = getSaturationUpwind('potential', 0, g, vT, T, mob, upw);

            nPh = numel(g);
            [upCells_v, upCells_G] = deal(repmat({cL}, 1, nPh));
            [x_v, s_v, x_G, s_G] = deal(cell(1, nPh));
            for phNo = 1:nPh
                
                upCells_v{phNo}(flag_v(:,phNo)) = cR(flag_v(:,phNo));
                x_v{phNo} = disc.transformCoords(x, upCells_v{phNo}, false, true);
                s_v{phNo} = disc.evaluateSaturation(x_v{phNo}, upCells_v{phNo}, sdof, state);
                
                upCells_G{phNo}(flag_G(:,phNo)) = cR(flag_G(:,phNo));
                x_G{phNo} = disc.transformCoords(x, upCells_G{phNo}, false, true);
                s_G{phNo} = disc.evaluateSaturation(x_G{phNo}, upCells_G{phNo}, sdof, state);
                
            end
            
            f_v = f(s_v{1}, 1-s_v{2}, upCells_v{1}, upCells_v{2});
            f_G = f(s_G{1}, 1-s_G{2}, upCells_G{1}, upCells_G{2});
            
            [x_c, ~, ~] = disc.transformCoords(x, cR);

            I = integrand(sdof, sdof, 1, 1, 1).*0;
            
            cellNo = G.mappings.cellMap.old2new(cR);
            f2c = sparse(cellNo, faceMap(faceNo), 1);
            
            cells = G.mappings.cellMap.old2new(bc.cell(:,2));
            
            sgn = 1 - 2*(bc.value>0);
            W = W.*sgn;
            
            for dofNo = 1:nDofMax
                
                keepCells = nDof(cells) >= dofNo;
                
                if any(keepCells)
                    
                    i  = W*integrand(f_v, s_v{1}, cellNo, faceNo, psi{dofNo}(x_c));
                    i  = f2c*i;
                    
                    ix = disc.getDofIx(state, dofNo, cells(keepCells)');
                    
                    I(ix) = i(cells(keepCells));
                    
                elseif numel(cells) == disc.G.cells.num
                    
                    warning('No cells with %d dofs', dofNo);
                    
                end
                
            end
            
            %

%             s = sL
%             
%             
%             sgn = 1 - 2*(G.faces.neighbors(faces, 1) ~= cells);
%             W = W.*sgn;
%             numFlux = nnz(isFlux);
%             
%             
%             
%             all2BC = nan(G.faces.num,1);
%             all2BC(faces) = 1:numFlux;
%             locFaceNo = all2BC(faceNo);
%             
%             sL = bc.sat(isFlux,1);
%             sL = sL(locFaceNo);
%             
%             isInj = bc.value(isFlux) > 0;
%             s = sL.*isInj(locFaceNo) + sR.*(~isInj(locFaceNo));
%             
%             cellNo = cR;
%             f = f(s, 1-s, cellNo, cellNo);
%             
%             I = integrand(sdof, sdof, 1, 1, 1).*0;
%             
%             for dofNo = 1:nDofMax
%                 
%                 keepCells = nDof(cells) >= dofNo;
%                 
%                 if any(keepCells)
%                     
%                     ix = disc.getDofIx(state, dofNo, cells(keepCells)');
%                     i  = W*integrand(f, s, cellNo, faceNo, psi{dofNo}(xR));
%                     I(ix) = i(keepCells);
%                     
%                 end
%                 
%             end
%             
            I = disc.trimValues(I);
            
        end
        
        function [flag_v, flag_G, upCells_v, upCells_G, s_v, s_G] = getUpstreamCell(disc, faces, x, T, vT, g, mob, sdof, state)
            
            G = disc.G;
            all2int = zeros(G.faces.num,1);
            all2int(disc.internalConn) = 1:nnz(disc.internalConn);
            ix = all2int(faces);
            
            cL = disc.N(ix,1);
            cR = disc.N(ix,2);
            T = T(ix);
            vT = vT(faces);
%             G = cellfun(@(g) g(ix), G, 'unif', false);
            g = cellfun(@(g) g(faces), g, 'unif', false);
            
            [xL, ~, ~] = disc.transformCoords(x, cL);
            [xR, ~, ~] = disc.transformCoords(x, cR);
            
%             sdof = double(sdof);
%             sL = disc.evaluateSaturation(xL, cL, sdof, state);
%             sR = disc.evaluateSaturation(xR, cR, sdof, state); 
            s = disc.evaluateSaturation([xL; xR], [cL; cR], sdof, state);
            sL = s(1:numel(cL));
            sR = s(numel(cL)+1:end);
            
            mob{1} = mob{1}([sL; sR], [cL; cR]);
            mob{2} = mob{2}(1-[sL; sR], [cL; cR]);
            
            N = [1:numel(ix); numel(ix)+1:2*numel(ix)]';
            
            upw = @(flag, x)faceUpstr(flag, x, N, [size(N,1), max(max(N))]);
            [flag_v, flag_G] = getSaturationUpwind('potential', 0, g, vT, T, mob, upw);
            
            nPh = numel(g);
            [upCells_v, upCells_G] = deal(repmat({cR}, 1, nPh));
            [s_v, s_G] = deal(cell(1, nPh));
            [s_v{:},s_G{:}] = deal(sR);
            for phNo = 1:nPh
                
                upCells_v{phNo}(flag_v(:,phNo)) = cL(flag_v(:,phNo));
                s_v{phNo}(flag_v(:,phNo)) = sL(flag_v(:,phNo));
%                 x_v{phNo} = disc.transformCoords(x, upCells_v{phNo});
%                 s_v{phNo} = disc.evaluateSaturation(x_v{phNo}, upCells_v{phNo}, sdof, state);
%                 
                upCells_G{phNo}(flag_G(:,phNo)) = cL(flag_G(:,phNo));
                s_G{phNo}(flag_G(:,phNo)) = sL(flag_G(:,phNo));
%                 x_G{phNo} = disc.transformCoords(x, upCells_G{phNo});
%                 s_G{phNo} = disc.evaluateSaturation(x_G{phNo}, upCells_G{phNo}, sdof, state);
                
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [W, x, cellNo, faceNo] = getCubature(disc, elements, type)
            
            if size(elements,1) == 1, elements = elements'; end
            
            G = disc.G;
            useMap = false;
            if isfield(G, 'mappings')
                useMap = true;
                maps = G.mappings;
                G = disc.G.parent;
                
                switch type 
                    case {'volume', 'surface'}
                        elements = maps.cellMap.new2old(elements);
                    case 'face'
                        elements = maps.faceMap.new2old(elements);
                end
                
            end
            
            switch type
                case 'volume'
                    
                    cubature = disc.volumeCubature;
                    ix = mcolon(cubature.parentPos(elements), cubature.parentPos(elements+1)-1);
%                     ixf = ix;
                    nq = diff(cubature.parentPos);
                    nq = nq(elements);
                    
                    cellNo = rldecode(elements, nq, 1);
                    faceNo = [];

                    sgn   = ones(numel(ix),1);
                    
                case 'face'
                    
                    cubature = disc.surfaceCubature;
                    nqf = diff(cubature.parentPos);
                    faceNo = rldecode(elements, nqf(elements), 1);
                    ix = mcolon(cubature.parentPos(elements), cubature.parentPos(elements+1)-1);
                    nq = nqf(elements);
                    cellNo = [];
                    sgn = 1;
                    
                case 'surface'
                    
                    cubature = disc.surfaceCubature;
                    faces = G.cells.faces(mcolon(G.cells.facePos(elements), G.cells.facePos(elements+1)-1));
                    ncf = accumarray(rldecode((1:G.cells.num)', diff(G.cells.facePos), 1), disc.internalConnParent(G.cells.faces(:,1)));
                    faces = faces(disc.internalConnParent(faces));
                    if size(faces,1) == 1, faces = faces'; end
                    ix = mcolon(cubature.parentPos(faces), cubature.parentPos(faces+1)-1);
                    
                    nqf = diff(cubature.parentPos);
                    nq = accumarray(rldecode((1:numel(elements))', ncf(elements), 1), nqf(faces));
                    if isempty(nq), nq = 0; end
                    
                    cellNo = rldecode(elements, nq);
                    faceNo = rldecode(faces, nqf(faces), 1);
                    sgn   = 1 - 2*(G.faces.neighbors(faceNo,1) ~= cellNo);
                    
%                     ixf = nan(numel(ix),1);
%                     ixf(disc.internalConn(faceNo)) = 1:nnz(disc.internalConn(faceNo));
                    
            end
            
            if useMap
                cellNo = maps.cellMap.old2new(cellNo);
                faceNo = maps.faceMap.old2new(faceNo);
            end
            
            x = cubature.points(ix, :);
            w = cubature.weights(ix).*sgn;
            [ii, jj] = blockDiagIndex(ones(numel(elements),1), nq);
            W = sparse(ii, jj, w);
%             W = cubature.W(cells, ixf);
            
        end
        
        %-----------------------------------------------------------------%
        function [smin, smax] = getMinMaxSaturation(disc, state)
            
            useMap = isfield(disc.G, 'mappings');
            if useMap
                maps = disc.G.mappings.cellMap;
                cells = find(maps.keep);
                G = disc.G.parent;
                faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1),1);
                s     = @(x, c) disc.evaluateSaturation(x, maps.old2new(c), state.sdof(:,1), state);                
            else
                G = disc.G;
                cells = (1:G.cells.num)';
                faces = G.cells.faces(:,1);
                s     = @(x, c) disc.evaluateSaturation(x, c, state.sdof(:,1), state);
            end
            
%             faces = G.cells.faces(:,1);
            nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
            
            nfn = diff(G.faces.nodePos);
            ncf = diff(G.cells.facePos);
            ncn = accumarray(rldecode((1:numel(cells))', ncf(cells), 1), nfn(faces));
            
            x = G.nodes.coords(nodes,:);
            
            cellNo = rldecode(cells, ncn, 1);
            x = disc.transformCoords(x, cellNo, false, true);
            
            s = s(x, cellNo);
            
            jj = rldecode((1:numel(cells))', ncn, 1);
            s = sparse(jj, (1:numel(s))', s);
            
            smax = full(max(s, [], 2));
            smin = full(min(s, [], 2));
            
        end
        
        %-----------------------------------------------------------------%
        function [jumpVal, faces, cells] = getInterfaceJumps(disc, state)
            
            G = disc.G;
            faces = find(disc.internalConn);
            cells = G.faces.neighbors(disc.internalConn,:);
            if isfield(G, 'mappings')
                order = G.mappings.cellMap.localOrder;
                isUpstr = any(order <= order(~G.cells.ghost)',2);
                keep = all(isUpstr(cells),2);
                faces = faces(keep);
            end

            cells = G.faces.neighbors(faces,:);                
            s     = @(x, c) disc.evaluateSaturation(x, c, state.sdof(:,1), state);                
            
            xF = G.faces.centroids(faces,:);
            
            cL  = cells(:,1);
            xFL = disc.transformCoords(xF, cL);
            cR  = cells(:,2);
            xFR = disc.transformCoords(xF, cR);
            
            jumpVal = abs(s(xFL, cL) - s(xFR, cR));
                        
            if 0
            ff = jumpVal > disc.jumpTolerance;
            faces = find(ff);
            cellsL = cL(faces);
            cellsR = cR(faces);
            for cNo = 1:numel(faces)
                
                figure(1); clf;
                
                subplot(2,2,1)
                disc.plotCellSaturation(state, cellsL(cNo));
                
                subplot(2,2,2)
                disc.plotCellSaturation(state, cellsR(cNo));
                
                subplot(2,2,3)
                plotGrid(G);
                plotGrid(G, cellsL(cNo), 'facec', 'r')
                
                subplot(2,2,4)
                plotGrid(G);
                plotGrid(G, cellsR(cNo), 'facec', 'r')
                
            end
            end
            
        end
        
        %-----------------------------------------------------------------%
        function state = limiter(disc, state)
            
            nDofMax = disc.basis.nDof;
            G       = disc.G;

            [jump, over, under] = deal(false(G.cells.num,1));
            
            if disc.jumpTolerance < Inf
                % Cells with interface jumps larger than threshold
                [jumpVal, ~, cells] = disc.getInterfaceJumps(state);
                
                j = accumarray(cells(:), repmat(jumpVal,2,1) > disc.jumpTolerance) > 0;
                jump = false(G.cells.num,1);
                jump(cells(:)) = j(cells(:));
                
                jump(state.degree == 0) = false;
                
            end
            
            if disc.outTolerance < Inf
                % Cells maximum/minimum value outside [-tol, 1+tol]
                [smin, smax] = disc.getMinMaxSaturation(state);
                over  = smax > 1 + disc.outTolerance;
                under = smin < 0 - disc.outTolerance;
            end
            
            bad = jump | over | under;
            state.outside = over | under;
            state.jump    = jump;
            
            if isfield(G, 'mappings')
                bad(G.cells.ghost) = false;
            end
                
            
%             state0 = state;
%             state.degree = repmat(disc.degree, G.cells.num, 1);
%             state = disc.mapDofs(state, state0);
%             state = disc.updateDisc(state);
            
            sdof = state.sdof;
           
            if any(bad)
                
                switch disc.limiterType
                    
                    case 'kill'
                        % Simple "limiter" that reduces to zero-degree for
                        % all bad cells

%                         state.degree(bad) = 0;
                        ix = disc.getDofIx(state, 1, bad);
%                         sdof(ix,:) = min(max(sdof(ix,:), 0), 1);
                        sdof(ix,:) = min(max(state.s(bad,:), 0), 1);
                        sdof(ix,:) = sdof(ix,:)./sum(sdof(ix,:),2);
                        state.s(bad,:) = sdof(ix,:);
                        ix = disc.getDofIx(state, 2:nDofMax, bad);
                        sdof(ix,:) = [];
                        state.sdof = sdof;
                        state.degree(bad) = 0;
                        
                    case 'adjust'
                        % Reduce to dG(1) and adjust slope
                        
%                         meanOutside = state.s(:,1) > 1 | state.s(:,1) < 0;
%                 
%                         ix = disc.getDofIx(state, 1, meanOutside);
%                         sdof(ix,:) = min(max(sdof(ix,:), 0), 1);
%                         sdof(ix,:) = sdof(ix,:)./sum(sdof(ix,:),2);
%                         state.degree(meanOutside) = 0;
%                         bad(meanOutside) = false;

                        for dNo = 1:G.griddim
                    
                            ix0 = disc.getDofIx(state, 1, under);
                            ix  = disc.getDofIx(state, dNo+1, under);
                            sdof(ix) = sign(sdof(ix)).*sdof(ix0);
%                             sdof(ix) = sign(sdof(ix)).*min(abs(sdof(ix0)), abs(sdof(ix)));

                            ix0 = disc.getDofIx(state, 1, over);
                            ix  = disc.getDofIx(state, dNo+1, over);
                            sdof(ix) = sign(sdof(ix)).*(1 - sdof(ix0));

                        end
                        state.degree(bad) = 1;
                        
                        ix = disc.getDofIx(state, (1+G.griddim+1):nDofMax, bad);
                        sdof(ix,:) = [];
%                         ix = disc.getDofIx(state, G.gr
                        state.sdof = sdof;
                        
                        
%                         both = over & under;
%                         
%                         
%                         
%                         %                 disc.dofPos = disc.updateDofPos();
%                 
%                         sdof = sdof(:,1);
%                         cells_o = find(over);
%                         cells_u = find(under);

                    case 'tvb'
                
                        error('Limiter type not implemented yet...');

                        
                end
                
%                 state = disc.getCellSaturation(state);
%                 ix = disc.getDofIx(state, 2:nDofMax, meanOutside | bad);
%                 sdof(ix,:) = [];
%                 state.sdof = s
                 
%                 
%                 % Reduce to first-order
%                 ix = disc.getDofIx(state, (1 + G.griddim + 1):nDofMax, cells);
%                 sdof(ix) = 0;
%                 
%                 % Reduce linear dofs so that saturaion is within [0,1]
% %                 ix0 = disc.getDofIx(1                , cells);
%                 for dNo = 1:G.griddim
%                     
%                     ix0 = disc.getDofIx(state, 1, cells_u);
%                     ix  = disc.getDofIx(state, dNo+1, cells_u);
%                     sdof(ix) = sign(sdof(ix)).*min(abs(sdof(ix0)), abs(sdof(ix)));
%                     
%                     ix0 = disc.getDofIx(state, 1, cells_o);
%                     ix  = disc.getDofIx(state, dNo+1, cells_o);
%                     sdof(ix) = sign(sdof(ix)).*(1 - sdof(ix0));
%                     
%                 end
%                 
%                 state.sdof(:,1) = sdof;
%                 state.sdof(:,2) = -sdof;
%                 
%                 ix0 = disc.getDofIx(state, 1, (1:G.cells.num)');
%                 state.sdof(ix0,2) = 1 - state.sdof(ix0,1);
%                 
%                 [smin, smax] = disc.getMinMaxSaturation(state);
%                 
%                 s = disc.getCellSaturation(state);
                
            end
            
        end
        
        %-----------------------------------------------------------------%
        function v = faceFlux2cellVelocity(disc, model, faceFlux)
            
            
            
        end
        
        %-----------------------------------------------------------------%
        function plotCellSaturation(disc, state, cellNo)
            
            G = disc.G;
            
            faces = G.cells.faces(G.cells.facePos(cellNo):G.cells.facePos(cellNo+1)-1);
            nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
            nodes = reshape(nodes, 2, [])';
            
            swap = G.faces.neighbors(faces,1) ~= cellNo;

            nodes(swap,:) = nodes(swap, [2,1]); nodes = nodes(:,1);
            
            x = G.nodes.coords(nodes,:);
            
            x = disc.transformCoords(x, cellNo);
            
            
            if disc.dim == 1
                
                n = 100;
                xx = linspace(-1,1,n)';
                xk = xx;
                
            elseif disc.dim == 2
                
                n = 10; 
                xx = linspace(-1, 1, n);
                [xx, yy] = ndgrid(xx);
                xx = [xx(:), yy(:)];
                
                [in, on] = inpolygon(xx(:,1), xx(:,2), x(:,1), x(:,2));
                keep = in;
                xk = xx(keep,:);
                
            elseif disc.dim == 3
                
            end
                
            cellNo = repmat(cellNo, size(xk,1), 1);
            s = disc.evaluateSaturation(xk, cellNo, state.sdof, state);
            
            if disc.dim > 1
                s = scatteredInterpolant(xk, s);
                s = reshape(s(xx(:,1), xx(:,2)), n,n)';
                surf(s);
            else
                plot(xx, s);
            end
            
        end
        
    end
        
end