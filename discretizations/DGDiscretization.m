classdef DGDiscretization < SpatialDiscretization
    
    properties

        degree              % Degree of discretization, dG(degree)
        degree0
        basis   = 'legendre' % Type of basis functions. Standard is tensor 
                            % products of Legendre polynomials.
        dim                 % Dimension of disc to facilitate e.g. 2D 
                            % simulations on horizontal slice of 3D
                            % reservoir
        
        useMomentFitting    % Bool to tell the method to use experimental 
                            % unstructured cubature class
        cellCubature        % Cubature for cell integrals
        faceCubature        % Cubature for face integrals
        
        jumpTolerance         % Tolerance for sat jumps across interfaces
        jumpLimiter
        outTolerance          % Tolerance for sat outside [0,1]
        outLimiter
        meanTolerance         % Tolerance of mean sat outside [0,1]
        limitAfterNewtonStep  % Use limiter after each Newton iteration
        limitAfterConvergence % Use limiter after convergence to reduce
                              % jumps and ensure values inside [0,1]
        
        plotLimiterProgress % 1d plot of result before and after limiter
        
        velocityInterp      % Function for mapping face fluxes to cell
                            % velocity/ies
        upwindType          % Type of upwind calculation
        
        internalConnParent  % If we only solve on subset of full grid, we
                            % must keep tract of internal connections in 
                            % the full grid.
                            
        limiter
                            
        nDof
        dofPos
        sample
        
        interp_setup
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function disc = DGDiscretization(G, varargin)
            
%             disc = disc@WENODiscretization(model, model.G.griddim, 'includeBoundary', true);
            disc = disc@SpatialDiscretization(G, varargin{:});
            
            G = disc.G;
            
            % Standard dG properties
            disc.degree = 1;
%             disc.degree0 = [];
            disc.basis  = 'legendre';
            disc.dim    = G.griddim;
            
            % Cubature
            disc.useMomentFitting = false;
            
            % Specifics for reordering
            disc.internalConnParent  = disc.internalConn;
            disc.G.parent            = G;
            
            [disc, basisArgs] = merge_options(disc, varargin{:});
            
            % Replace basis string by dG basis object
            if ~isfield(disc.basis, 'psi')
                disc.basis  = dgBasis(disc.dim, disc.degree, disc.basis, basisArgs{:});
            end
            disc.degree = disc.basis.degree;
%             disc.nDof   = disc.getnDof();
%             [~, disc] = disc.updateDofPos(disc);
            
            % Set up velocity interpolation
            disc.velocityInterp = velocityInterpolation(G, 'mimetic');
            disc.upwindType     = 'potential';
            
            % Create cubatures
            prescision = 2*max(disc.degree);
            isCoarse   = isfield(G, 'parent');
            if G.griddim == 2
                if isCoarse %&& ~disc.useUnstructCubature
                    volCub  = CoarseGrid2DCubature(G, prescision);
                else
                    if all(disc.degree == 0) || disc.useMomentFitting
                        volCub = MomentFitting2DCubature(G, prescision);
                    else
                        volCub = TriangleCubature(G, prescision);
                    end
                end
                surfCub = LineCubature(G, prescision);
            else
                if isCoarse 
                    if ~disc.useMomentFitting
                        volCub  = CoarseGrid3DCubature(G, prescision);
                        surfCub = TriangleCubature(G, prescision);
                    else
                        volCub  = MomentFitting3DCubature(G, prescision);
                        surfCub = TriangleCubature(G, prescision);
                    end
                else
                    if all(disc.degree == 0)
                        volCub  = MomentFitting3DCubature(G, prescision);
                        surfCub = MomentFitting2DCubature(G, prescision);
                    elseif disc.useMomentFitting
                        volCub  = MomentFitting3DCubature(G, prescision);
                        surfCub = TriangleCubature(G, prescision);
                    else
                        volCub  = TetrahedronCubature(G, prescision);
                        surfCub = TriangleCubature(G, prescision);
                    end
                end
            end
            
            disc.cellCubature  = volCub;
            disc.faceCubature = surfCub;
            
        end
        
        %-----------------------------------------------------------------%
        function [state, disc] = updateDofPos(disc, state)
            % Update dosfPos (position of dofs in state.sdof) based on
            % changes in state.degree, or create dofPos vector if it does
            % not exist. Dofs for cell i are found in
            %
            %   state.sdof(dofPos(:,i),:),
            %
            % Zeros are included to easily map dofs from one timestep to
            % the next.

            dp = reshape((1:disc.G.cells.num*disc.basis.nDof)', disc.basis.nDof, []);
            if isfield(state, 'nDof')
                nd = state.nDof;
            else
                nd = disc.getnDof(state);
            end
            subt = cumsum([0; disc.basis.nDof - nd(1:end-1)]);
            [ii, jj, v] = find(dp);

            if size(ii,1) == 1, ii = ii'; end
            if size(jj,1) == 1, jj = jj'; end
            if size(v ,1) == 1, v  =  v'; end

            cnDof = cumsum(nd);

            v = v - subt(jj);
            v(v > cnDof(jj)) = 0;
            dp = full(sparse(ii, jj, v));
            
            state.nDof   = nd;
            state.dofPos = dp;
            disc.nDof    = nd;
            disc.dofPos  = dp;
            
        end
        
        %-----------------------------------------------------------------%
        function ix = getDofIx(disc, state, dofNo, cells, includezero)
            % Get position of dofs in state.sdof for a given cell
            %
            % PARAMETERS:
            %   state       - State with field sdof
            %   dofNo       - Dof number we want the position of. Empty
            %                 dofNo returns position of all dofs for cells
            %   cells       - Cells we want the dof position for. If empty,
            %                 positions of dofNo are returned for all cells.
            %   includeZero - Boolean indicating of we should include zeros
            %                 or not (see updateDofPos)
            %
            % RETURNS:
            %   ix - Indices into states.sdof. Dof number dofNo for cells
            %        are found in state.sdof(ix,:);

            G = disc.G;
            if nargin < 3 || (numel(dofNo) == 1 && dofNo == Inf)
                % dofNo not given, return ix for all dofs
                dofNo = 1:disc.basis.nDof;
                cells = 1:G.cells.num;
            elseif nargin < 4 || (numel(cells) == 1 && cells == Inf)
                % Cells not given, return ix for all cells
                cells = 1:G.cells.num;
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
        function [xhat, translation, scaling] = transformCoords(disc, x, cells, inverse, useParent)
            % Transfor coordinates from physical to reference coordinates. 
            %
            % PARAMETERS:
            %   x         - Coordinates in physical space
            %   cells     - Cells we want reference coordinates for, cells(ix)
            %               are used when transforming x(ix,:)
            %   inverse   - Boolean indicatiing if we are mapping 
            %               to (inverse = false) or from (inverse = true)
            %               reference coordiantes. Default = false.
            %   useParent - Boolean indicating if we are working on the
            %               full grid (G.parent) or a subgrid.
            %
            % RETURNS:
            %   xhat        - Transformed coordinates
            %   translation - Translation applied to x
            %   scaling     - Scaling applied to x
            
            G = disc.G;
            
            if nargin < 4, inverse   = false; end
            if nargin < 5, useParent = false; end
            
            if isfield(G, 'mappings') && useParent
                G = G.parent;
            end
            
            [xhat, translation, scaling] = disc.cellCubature.transformCoords(x, cells, inverse);
            
            if 0
            % Coordinates are centered in cell center
            if any(strcmpi(G.type, 'generateCoarseGrid'))
                translation = -G.cells.centers(cells,:);
            else
                translation = -G.cells.centroids(cells,:);
            end
            if isfield(G.cells, 'dx')
                % Scaling found from dimensions of minimum bounding box
                % aligned with coordinate axes that contains the cell
                scaling = 1./(G.cells.dx(cells,:)/2);
            else
                % If it G.cells.dx is not computed, we use approximation
                dx = G.cells.volumes(cells).^(1/G.griddim);
                scaling = repmat(1./(dx/2), 1, disc.dim);
            end
            
            if ~inverse
                xhat = (x + translation).*scaling;
                xhat = xhat(:, 1:disc.dim);
                scaling     = scaling(:, 1:disc.dim);
                translation = translation(:, 1:disc.dim);
%                 assert(all(all(abs(xhat)<=1)))
            else
                xhat = x./scaling - translation;
            end
               
            end
        end
        
        %-----------------------------------------------------------------%
        function nDof = getnDof(disc, state)
            % Get number of dofs for each cell from degree.
            if nargin == 1
                deg = repmat(disc.degree, disc.G.cells.num, 1);
            else
                deg = state.degree;
            end
            if disc.degree < 0
                nDof = 0;
            else
                nDof = factorial(deg + disc.dim)...
                       ./(factorial(disc.dim).*factorial(deg));
            end
            
        end
        
        %-----------------------------------------------------------------%
        function state = mapDofs(disc, state, state0, name)
            % Map dofs from state0 to state, typically from one timestep to
            % the next, when we start with maximum number of dofs in all
            % cells.
            
            % Update dofPos
            state = disc.updateDofPos(state);
            
            if nargin == 3
                name = 'sdof';
            else
                name = [name, 'dof'];
            end
            if all(state.nDof == state0.nDof)
                state.(name) = state0.(name);
            else
                dof = zeros(sum(state.nDof), size(state0.(name),2));
                for dofNo = 1:disc.basis.nDof
                    % We may be solving only a subset of the gridcells, so
                    % we include zeros in dofIx to keep track of where old
                    % dofs maps to new ones.
                    ix  = disc.getDofIx(state , dofNo, (1:disc.G.cells.num)', true);
                    ix0 = disc.getDofIx(state0, dofNo, (1:disc.G.cells.num)', true);
                    dof(ix(ix0 > 0 & ix > 0),:) = state0.(name)(ix0(ix0 > 0 & ix > 0),:);
                end
                state.(name) = dof;
            end

        end
        
        %-----------------------------------------------------------------%
        function p = evaluateDGVariable(disc, x, cells, state, pdof, psi)
            
            nDof    = state.nDof; %#ok
            nDofMax = disc.basis.nDof;
            
            if isempty(pdof)
                return
            else
                if nargin < 6 || isempty(psi)
                    x    = disc.transformCoords(x, cells);
                    psi  = disc.basis.psi;
                    getx = @(keep) x(keep,:);
                else
                    getx = @(keep) keep;
                end
                p = pdof(cells)*0;
                for dofNo = 1:nDofMax
                    keep = nDof(cells) >= dofNo; %#ok
                    ix = disc.getDofIx(state, dofNo, cells(keep));
                    if all(keep)
                        p = p + pdof(ix,:).*psi{dofNo}(getx(keep));
                    else
                        p(keep, :) = p(keep, :) + pdof(ix,:).*psi{dofNo}(getx(keep));
                    end
                end
            end

        end
        
        %-----------------------------------------------------------------%
        function p = evaluateProp(disc, state, dof, type, elements)
            psi = [];
            if nargin < 5
                elements = Inf;
            end
            switch type
                case 'cell'
                    [~ , x, cNo] = disc.getCubature(elements, 'cell');
                    if isfield(state, 'psi_c') && all(elements < Inf)
                        psi = state.psi_c;
                    end
                case 'face'
                    [~ , x, ~, fNo] = disc.getCubature(elements, 'face');
                    x   = repmat(x, 2, 1);
                    N   = disc.G.faces.neighbors;
                    cNo = [N(fNo,1); N(fNo,2)];
                    keep = cNo ~= 0;
                    cNo = cNo(keep);
                    x = x(keep,:);
                    if isfield(state, 'psi_f') && all(elements < Inf)
                        psi = state.psi_f;
                    end
            end
            if ~isempty(psi) && ~iscell(psi)
                psi = mat2cell(psi, size(psi,1), ones(1,disc.basis.nDof));
            end
            if iscell(dof)
                p = cell(1,numel(dof));
                for i = 1:numel(dof)
                    p{i} = disc.evaluateDGVariable(x, cNo, state, dof{i}, psi);
                end
            else
                p = disc.evaluateDGVariable(x, cNo, state, dof, psi);
            end
        end
        
        function fill = getFillSat(disc, state)
            fill = zeros(sum(state.nDof),1);
            ix   = disc.getDofIx(state, 1);
            fill(ix) = 1;
        end
        
        %-----------------------------------------------------------------%
        function varargout = evaluateDGVariables(disc, x, cells, state, varargin)
            varargout = cellfun(@(v) disc.evaluateDGVariable(x, cells, state, v), ...
                                               varargin, 'UniformOutput', false);
        end
        
        %-----------------------------------------------------------------%
        function varargout = getCellMean(disc, state, varargin)
            % Get average cell value from dofs

            % Get cubature for all cells
            W = disc.getCubature(Inf, 'cell');
            val = cell(numel(varargin),1);
            
            for i = 1:nargin-2
                dof = varargin{i};
                v = disc.evaluateProp(state, dof, 'cell');
                if iscell(dof)
                    for j = 1:numel(dof)
                        v{j} = W*v{j};
                    end
                else
                    v = W*v;
                end
                val{i} = v;
            end
            varargout = val;
            
        end
        
        %-----------------------------------------------------------------%
        function dotProduct = dot(disc,u,v)
            dotProduct = sum(u.*v, 2);
        end
        
        %-----------------------------------------------------------------%
        function ip = inner(disc, u, v, differential, elements, bc)
            
            if nargin < 5
                elements = Inf;
            end
                
            switch differential
                case 'dV'
                    ip = disc.cellInt(u, v, elements);
                case 'dS'
                    ip = disc.faceInt(u, v, elements);
                case 'dSbc'
                    ip = disc.faceFluxIntBC2(u, v, bc);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function I = cellInt(disc, u, v, cells)
            % Integrate integrand over cells
            %
            % PARAMETERS:
            %   model    - Model, which contains information on how the
            %              integrand looks like
            %   fun      - Integrand function handle
            %   cells    - Cells over which we will integrate fun
            %   state, state0 - States from current and prev timestep,
            %              to be used for dofPos
            %   varargin - Variables passed to model for integrand
            %              evalauation. varargin{1} MUST be an AD object of
            %              dofs for all cells in the grid.
            %
            % RETURNS:
            %   I - Integrals int(fun*psi{dofNo}) for dofNo = 1:nDof over
            %       all cells
            
            nDofMax  = numel(v);     % Maximum number of dofs
            % Empty cells means all cells in grid
            if isempty(cells)
                cells = (1:disc.G.cells.num)';
            end
            % Get cubature for all cells, transform coordinates to ref space
            [W, x, cellNo, ~] = disc.getCubature(cells, 'cell');
            if isinf(cells)
                cells = (1:disc.G.cells.num)';
            end
            [x, ~, scaling]   = disc.transformCoords(x, cellNo);
            % Evaluate integrals
%             I = dof*0;
            I = disc.sample*0;
            for dofNo = 1:nDofMax
                keepCells = disc.nDof(cells) >= dofNo;
                if any(keepCells)
                    ix = disc.getDofIx(disc, dofNo, cells(keepCells));
                    if isa(u, 'SpatialVector')
                        i = W*disc.dot(u,v{dofNo}(x).*scaling);
                    else
                        i = W*(u.*v{dofNo}(x));
                    end
                    I(ix) = i(keepCells);
                elseif numel(cells) == disc.G.cells.num
                    warning('No cells with %d dofs', dofNo);
                end
            end
            I = disc.trimValues(I);
        end
        
        %-----------------------------------------------------------------%
        function I = faceInt(disc, u, v, faces)
            % Integrate integrand over all internal faces of each cell in
            % cells
            %
            % PARAMETERS:
            %   model    - Model, which contains information on how the
            %              integrand looks like
            %   fun      - Integrand function handle
            %   cells    - Cells over which we will integrate fun
            %   state    - State to be used for dofPos
            %   varargin - Variables passed to model for integrand
            %              evalauation. varargin{1} MUST be an AD object of
            %              dofs for all cells in the grid.
            %
            % RETURNS:
            %   I - Integrals int(fun*psi{dofNo}) for dofNo = 1:nDof over
            %       all cell surfaces of all cells
            
            nDofMax = disc.basis.nDof; % maximum number of dofs
            % Get cubature
            [W, x, ~, faceNo] = disc.getCubature(faces, 'face');
            if isinf(faces)
                % faces = inf means all internal faces in the grid
                faces = find(disc.internalConn)';
            end
            N = disc.G.faces.neighbors;
            if isempty(faceNo)
                I = 0;
                return
            end
            I = disc.sample*0;
            facesSide = faces;
            xSide     = x;
            uSide     = u;
            WSide     = W;
            for side = 1:2
                cells  = N(faces , side);
                cellNo = N(faceNo, side);
                ix_c   = cells > 0;
                if ~all(ix_c)
                    cells     = cells(ix_c);
                    facesSide = faces(ix_c);
                    ix_cNo    = cellNo > 0;
                    cellNo    = cellNo(ix_cNo);
                    xSide     = x(ix_cNo,:);
                    uSide     = u(ix_cNo);
                    WSide     = W(ix_c,ix_cNo);
                end
                f2c = sparse(cells, (1:numel(facesSide))', 1, disc.G.cells.num, numel(facesSide));
                [xf, ~, ~] = disc.transformCoords(xSide, cellNo);
                % Evaluate integrals
                sgn = (-1).^(side-1);
                for dofNo = 1:nDofMax
                    keepCells = disc.nDof(cells) >= dofNo;
                    kc = f2c*keepCells > 0;
                    if any(keepCells)
                        ix = disc.getDofIx(disc, dofNo, kc');
                        i  = f2c*(sgn.*WSide*(uSide.*v{dofNo}(xf)));
                        I(ix) = I(ix) + i(kc);
                    elseif numel(facesSide) == nnz(disc.internalConn)
                        warning('No cells with %d dofs', dofNo);
                    end
                end
                I = disc.trimValues(I);
            end
        end
        
         %-----------------------------------------------------------------%
        function I = faceFluxIntBC2(disc, u, v, bc)
            % Integrate integrand over all faces where bcs are defined
            %
            % PARAMETERS:
            %   model    - Model, which contains information on how the
            %              integrand looks like
            %   fun      - Integrand function handle
            %   bc       - Boundary condition struct from schedule
            %   state    - State to be used for dofPos
            %   varargin - Variables passed to model for integrand
            %              evalauation. varargin{1} MUST be an AD object of
            %              dofs for all cells in the grid.
            %
            % RETURNS:
            %   I - All integrals int(fun*psi{dofNo}) for dofNo = 1:nDof
            %       over all bc faces for each cell
            
            G       = disc.G;          % Grid
            nDofMax = disc.basis.nDof; % Maximum number of dofs
            % Get faces and corresponding cells where BCs are defined
            faces = bc.face;
            cells = sum(G.faces.neighbors(faces,:),2);
            % Get cubature for each face, find corresponding cells, and
            % transform to reference coords
            [W, x, ~, faceNo] = disc.getCubature(faces, 'face');
            cellNo = sum(G.faces.neighbors(faceNo,:),2);
            [x_r, ~, ~] = disc.transformCoords(x, cellNo);
            % Ensure that we have the right sign for the integrals
            sgn = 1 - 2*(G.faces.neighbors(faces, 1) == 0);
            W = W.*sgn;
            % Mappings from global cell numbers to bc face/cell numbers
            globCell2BCcell = nan(G.cells.num,1);
            globCell2BCcell(cells) = 1:numel(cells);
            S = sparse(globCell2BCcell(cells), 1:numel(faces), 1);
            % Evaluate integrals
            I = disc.sample*0;
            for dofNo = 1:nDofMax
                keepCells = disc.nDof(cells) >= dofNo;
                if any(keepCells)
                    ix = disc.getDofIx(disc, dofNo, cells(keepCells)');
                    i  = W*(u.*v{dofNo}(x_r));
                    i  = S*i;
                    I(ix) = i(keepCells);
                end
            end
            %I = disc.trimValues(I);
            
        end
        
        %-----------------------------------------------------------------%
        function I = faceFluxIntBC(disc, fun, bc, state, sdof)
            % Integrate integrand over all faces where bcs are defined
            %
            % PARAMETERS:
            %   model    - Model, which contains information on how the
            %              integrand looks like
            %   fun      - Integrand function handle
            %   bc       - Boundary condition struct from schedule
            %   state    - State to be used for dofPos
            %   varargin - Variables passed to model for integrand
            %              evalauation. varargin{1} MUST be an AD object of
            %              dofs for all cells in the grid.
            %
            % RETURNS:
            %   I - All integrals int(fun*psi{dofNo}) for dofNo = 1:nDof
            %       over all bc faces for each cell
            
            G       = disc.G;          % Grid
            psi     = disc.basis.psi;  % Basis functions
            nDof    = state.nDof;      % Number of dofs per cell
            nDofMax = disc.basis.nDof; % Maximum number of dofs
            % Get faces and corresponding cells where BCs are defined
            faces = bc.face;
            cells = sum(G.faces.neighbors(faces,:),2);
            % Get cubature for each face, find corresponding cells, and
            % transform to reference coords
            [W, x, ~, faceNo] = disc.getCubature(faces, 'face');
            cellNo = sum(G.faces.neighbors(faceNo,:),2);
            [xR, ~, ~] = disc.transformCoords(x, cellNo);
            % Ensure that we have the right sign for the integrals
            sgn = 1 - 2*(G.faces.neighbors(faces, 1) == 0);
            W = W.*sgn;
            % Mappings from global cell numbers to bc face/cell numbers
            globCell2BCcell = nan(G.cells.num,1);
            globCell2BCcell(cells) = 1:numel(cells);
            S = sparse(globCell2BCcell(cells), 1:numel(faces), 1);
            % Evaluate integrals
%             I = getSampleAD(sdof)*0;
%             I = getSampleAD(sdof);
            I = sdof;
            for dofNo = 1:nDofMax
                keepCells = nDof(cells) >= dofNo;
                if any(keepCells)
                    ix = disc.getDofIx(state, dofNo, cells(keepCells)');
                    i  = W*fun(psi{dofNo}(xR));
                    i = S*i;
                    I(ix) = i(keepCells);
                end
            end
            I = disc.trimValues(I);
            
        end
        
        %-----------------------------------------------------------------%
        function [flag_v, flag_G, upCells_v, upCells_G, s_v, s_G] = getSaturationUpwind(disc, faces, x, T, vT, state, g, mob, sdof, rdof, cdof)
            % Explicit calculation of upstream cells. See getSaturationUpwindDG
            [flag_v, flag_G, upCells_v, upCells_G, s_v, s_G] ...
                = getSaturationUpwindDG(disc, faces, x, T, vT, state, g, mob, sdof, rdof, cdof);
        end
        
        %-----------------------------------------------------------------%
        function flag = multiphaseUpwindIndices(disc, G, v, T, mob, upstr)
            % Explicit calculation of upstream cells. See getSaturationUpwindDG
            
            [~, ~, ~, faces] = disc.getCubature(find(disc.internalConn), 'face');
            nf = numel(faces);
            % Mapping from all faces to internal connections
            all2int = zeros(disc.G.faces.num,1);
            all2int(disc.internalConn) = 1:nnz(disc.internalConn);
            ix = all2int(faces);
            
            T  = T(ix);
            v  = v(ix);
            
            % Make fake faceUpstr function
            upw = @(flag, x)faceUpstr(flag, x, [1:nf; nf+1:2*nf]', [nf, 2*nf]);

            flag = multiphaseUpwindIndices(G, v, T, mob, upw);
        end
        
        % --------------------------------------------------------------------%
        function gdxyz = getGravityGradient(disc, model)
            
            x = model.G.cells.centroids;
            [~, ~, ~, faces] = disc.getCubature(find(model.operators.internalConn), 'face');
            cells = [model.G.faces.neighbors(faces,1); model.G.faces.neighbors(faces,2)];
            gdxyz = model.operators.Grad(model.G.cells.centroids) * g';
        end
        
        %-----------------------------------------------------------------%
        function [W, x, cellNo, faceNo] = getCubature(disc, elements, type, varargin)
            % Get cubature for elements. Wrapper for cubature class
            % function getCubature, with mapping of elements before and
            % after in case we are solving on a subgrid
            
            useMap = isfield(disc.G, 'mappings');
            if useMap
                % Map elements to old numbering
                maps = disc.G.mappings; 
                switch type 
                    case {'cell', 'surface'}
                        elements = maps.cellMap.new2old(elements);
                    case 'face'
                        elements = maps.faceMap.new2old(elements);
                end
            end
            
            % Get correct cubature type
            switch type 
                case 'cell'
                    cubature = disc.cellCubature; 
                case {'surface', 'face'}
                    cubature = disc.faceCubature;
            end
            
            opt = struct('excludeBoundary', true                   , ...
                         'internalConn'   , disc.internalConnParent, ...
                         'outwardNormal'  , true                   );
            opt = merge_options(opt, varargin{:});
            
            % Get cubature from cubature class
            [W, x, ~, cellNo, faceNo] = cubature.getCubature(elements, type, ...
                                 'excludeBoundary', opt.excludeBoundary    , ...
                                 'internalConn'   , opt.internalConn       , ...
                                 'outwardNormal'  , opt.outwardNormal      );
                
            if useMap
                % Map elements back to new numbering
                cellNo = maps.cellMap.old2new(cellNo);
                faceNo = maps.faceMap.old2new(faceNo);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [sMin, sMax] = getMinMaxSaturation(disc, state)
            % Get maximum and minimum saturaiton for each cell
            
            G = disc.G;
            cells = (1:G.cells.num)';
            
            % Get all quadrature points for all cells
            [~, xSurf, cSurf, ~] = disc.getCubature(cells, 'surface', 'excludeBoundary', false);
            [~, xCell, cCell, ~] = disc.getCubature(cells, 'volume' );
            % Evaluate saturation at faces
            sSurf = disc.evaluateDGVariable(xSurf, cSurf, state, state.sdof(:,1));
            [~, nSurf] = rlencode(cSurf);
            [sMins, sMaxs] = getMinMax(sSurf, nSurf);
            % Evaluate saturation at cells
            sCell = disc.evaluateDGVariable(xCell, cCell, state, state.sdof(:,1));
            [~, nCell] = rlencode(cCell);
            % Find min/max saturation
            [sMinc, sMaxc] = getMinMax(sCell, nCell);
            sMin = min(sMins, sMinc);
            sMax = max(sMaxs, sMaxc);
            
        end
        
        %-----------------------------------------------------------------%
        function [vMin, vMax] = getMinMax(disc, state, dof)
            % Get maximum and minimum saturaiton for each cell
            
            G = disc.G;
            
            % Get all quadrature points for all cells
            [~, xF, ~, fF] = disc.getCubature((1:G.faces.num), 'face');
            xF = repmat(xF, 2, 1);
            cF = [G.faces.neighbors(fF,1); G.faces.neighbors(fF,2)];
            ix = cF == 0;
            cF(ix) = [];
            xF(ix,:) = [];
            [cF, ix] = sort(cF);
            xF = xF(ix,:);
            [~, xC, cC, ~] = disc.getCubature(Inf, 'cell');
            % Evaluate saturation at faces
            vF = disc.evaluateDGVariable(xF, cF, state, dof);
            [~, nF] = rlencode(cF);
            [vMinF, vMaxF] = getMinMax(vF, nF);
            % Evaluate saturation at cells
            vC = disc.evaluateDGVariable(xC, cC, state, dof);
            [~, nC] = rlencode(cC);
            % Find min/max saturation
            [vMinC, vMaxC] = getMinMax(vC, nC);
            vMin = min(vMinF, vMinC);
            vMax = max(vMaxF, vMaxC);
            
        end
        
        %-----------------------------------------------------------------%
        function [jumpVal, faces, cells] = getInterfaceJumps(disc, dof, state)
            % Get interface jumps for all internal connections in grid
            
            G     = disc.G;
            faces = find(disc.internalConn);
            cells = G.faces.neighbors(disc.internalConn,:);
            if isfield(G, 'mappings')
                % We assume we are using reordering, and thus only check
                % interface jumps for interfaces against already solved
                % cells
                order   = G.mappings.cellMap.localOrder;
                isUpstr = any(order <= order(~G.cells.ghost)',2);
                keep    = all(isUpstr(cells),2);
                faces   = faces(keep);
            end
            
            % Saturation function
            s = @(x, c) disc.evaluateDGVariable(x, c, state, dof);
            
            % Get reference coordinates
            xF    = G.faces.centroids(faces,:);            
            cells = G.faces.neighbors(faces,:);                
            cL    = cells(:,1);
            cR    = cells(:,2);
            
            % Find inteface jumps
            jumpVal = abs(s(xF, cL) - s(xF, cR));
            
        end
        
        %-----------------------------------------------------------------%
        function state = limiters(disc, model, state, state0, before)
            % Limiters applied after each Newton iteration if
            % disc.limitAfterNewtonStep = true and before = true, and after
            % convergence if disc.limitAfterConvergence = true and before =
            % false

            G = disc.G;
            check = true(G.cells.num,1);
            if isfield(G, 'mappings')
                check(G.cells.ghost) = false;
            end
                
            if before
                
                % Limiters to be applied after each Newton iteration
%                 if disc.meanTolerance < Inf
%                     % If the mean cell saturation is outside [0,1],
%                     % something is wrong, and we reduce to order 0
%                     s = state.s;
%                     meanOutside = any(s < 0 - disc.meanTolerance | ...
%                                       s > 1 + disc.meanTolerance, 2);
%                     bad = meanOutside & check;
%                     if any(bad)
%                         state = dgLimiter(disc, state, bad, 'kill', 'plot', disc.plotLimiterProgress);
%                     end 
%                 end
                
                if disc.degree > 0
                    if disc.outTolerance < Inf && 1
                        % If saturation is outside [0,1], we reduce to order 1
                        if 0
                            state = dgLimiter(disc, state, check, 's', 'scale', 'plot', disc.plotLimiterProgress);
                        else
                            [sMin, sMax] = disc.getMinMaxSaturation(state);
                            outside = sMin < 0 - disc.outTolerance | ...
                                      sMax > 1 + disc.outTolerance;
                            bad = outside & check;
                            if any(bad)
                                state = dgLimiter(disc, state, bad, 's', disc.outLimiter, 'plot', disc.plotLimiterProgress);
                            end
                        end
                    end

                    if disc.jumpTolerance < Inf && 1
                        % Cells with interface jumps larger than threshold
                        [jumpVal, ~, cells] = disc.getInterfaceJumps(state.sdof(:,1), state);
                        j = accumarray(cells(:), repmat(jumpVal,2,1) > disc.jumpTolerance) > 0;
                        jump = false(G.cells.num,1);
                        jump(cells(:))          = j(cells(:));
                        jump(state.degree == 0) = false;
%                         if size(jump,1) == 1; jump = jump'; end
                        bad = jump & check;
                        if any(bad)
                            state = dgLimiter(disc, state, bad, 's', disc.jumpLimiter, 'plot', disc.plotLimiterProgress);
                        end
                    end
                end
            else
                % Limiters to be applied after convergence
                if disc.degree > 0
                    % Scale solutions so that 0 <= s <= 1
%                     state = dgLimiter(disc, state, check, 's', 'scale', 'plot', disc.plotLimiterProgress);
                    if disc.jumpTolerance < Inf
                        % Limit saturation slope in cells with interface jumps
                        % larger than threshold
                        [jumpVal, ~, cells] = disc.getInterfaceJumps(state.sdof(:,1), state);
                        j = accumarray(cells(:), repmat(jumpVal,2,1) > disc.jumpTolerance) > 0;
                        jump = false(G.cells.num,1);
                        jump(cells(:))          = j(cells(:));
                        jump(state.degree == 0) = false;
                        bad = jump & check & state.degree > 0;
                        if any(bad)
                            state = dgLimiter(disc, state, bad, 's', 'tvb', 'plot', disc.plotLimiterProgress);
                        end
                    end
                end
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
        
    end
        
end