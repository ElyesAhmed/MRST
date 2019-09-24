classdef MomentFitting2DCubature < Cubature
    % Cubature based on moment-fitting for MRST grids
    
    properties
        
        reduce
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = MomentFitting2DCubature(G, prescision, internalConn, varargin)
            % Set up cubature
            
            % Basic properties handled by parent class
            cubature = cubature@Cubature(G, prescision, internalConn);
            cubature.reduce = true;
            cubature.dim       = 2;
            cubature = merge_options(cubature, varargin{:});
            % Make cubature points and weights
            [x, w, n] = cubature.makeCubature();
            % Assing properties
            cubature.points    = x;
            cubature.weights   = w;
            cubature.numPoints = n;
            % Construct cubature position vector
            cubature.pos = [0; cumsum(n)] + 1;
            
        end
           
        %-----------------------------------------------------------------%
        function [x, w, n] = makeCubature(cubature)
            
            % Dimension of cubature
            dim = 2;
            G   = cubature.G;
            % Basis functions used in moment-fitting
            basis = dgBasis(dim, cubature.prescision, 'legendre');
            psi   = basis.psi;
            nDof  = basis.nDof;
            
            if cubature.prescision <= 1
                % Precision is 1, use midpoint rule
                [x, ~, n] = getSquareCubaturePointsAndWeights(cubature.prescision);
                if G.griddim == 2
                    w = G.cells.volumes;
                    x = repmat(x, G.cells.num, 1);
                    n = repmat(n, G.cells.num, 1);
                    type = 'volume';
                else
                    w = G.faces.areas;
                    type = 'face';
                    x = repmat(x, G.faces.num, 1);
                    n = repmat(n, G.faces.num, 1);
                end
            else
                % If precision is higher than one, we must calculate
                % weights based on moment-fitting
                
                % The starting point is a quadrature for a reference square
                % with more than nDof points
                if 1
                    G1 = computeGeometry(cartGrid([1,1], [2,2]));
                    G1.nodes.coords = G1.nodes.coords - 1;
                    G1 = computeVEMGeometry(G1);
                    G1 = computeCellDimensions2(G1);
                    cubTri = TriangleCubature(G1, cubature.prescision, cubature.internalConn);
                    [~, x, ~, cellNo, ~] = cubTri.getCubature(1, 'volume');
                    x = cubTri.transformCoords(x, cellNo);
                    x = unique(x, 'rows');
                else
                    k = 0;
                    nS = 0;
                    while nS < 2*nDof
                        [x, ww, nS] = getSquareCubaturePointsAndWeights(cubature.prescision + k);
                        k = k+1;
                    end
                end
                % Cubature is either for faces or cells
                if G.griddim > cubature.dim
                    type = 'face';
                    elements = 1:G.faces.num;
                else
                    type = 'volume';
                    elements = 1:G.cells.num;
                end
                % We use known cubature to calculate the moments
                if isfield(G, 'parent')
                    knownCub = CoarseGrid2DCubature(G, cubature.prescision, cubature.internalConn);
                else
                    knownCub = TriangleCubature(G, cubature.prescision, cubature.internalConn);
                end
                [~, xq, wq, cellNo, faceNo] = knownCub.getCubature(elements, type);
                switch type
                    case 'face'
                        wq = wq./G.faces.areas(faceNo);
                    case 'volume'
                        wq = wq./G.cells.volumes(cellNo);
                end
                % Map cubature points to reference coordinates
                if G.griddim == 3
                    % Map to face reference coordinates
                    vec1  = G.faces.coordSys{1}(faceNo,:);
                    vec2  = G.faces.coordSys{2}(faceNo,:);
                    xq    = xq - G.faces.centroids(faceNo,:);
                    xq    = [sum(xq.*vec1,2), sum(xq.*vec2, 2)];
                    xq    = xq./(G.faces.dx(faceNo,:)/2);
                    count = faceNo;
                    num   = G.faces.num;
                else
                    % Map to cell reference coordiantes
                    xq    = cubature.transformCoords(xq, cellNo);
                    count = cellNo;
                    num   = G.cells.num;
                end
                % Moments
                M = cellfun(@(p) accumarray(count, wq.*p(xq)), psi, 'unif', false);
                % Compute right-hand side
                rhs = zeros(nDof, num);
                tol = eps(mean(G.cells.volumes));
                for dofNo = 1:nDof
                    m = M{dofNo};
                    m(abs(m) < tol) = 0;
                    rhs(dofNo, :) = m;
                end
                [x,w,n] = fitMoments3(x, basis, M, 'reduce', cubature.reduce);
            end
            
            % Map from reference to physical coordinates
            if strcmp(type, 'face')
                % Face coordinates
                faceNo = rldecode((1:G.faces.num)', n, 1);
                vec1   = G.faces.coordSys{1}(faceNo,:);
                vec2   = G.faces.coordSys{2}(faceNo,:);
                dx     = G.faces.dx;
                x = x.*dx(faceNo,:)/2;                
                x = x(:,1).*vec1 + x(:,2).*vec2;
                x = x + G.faces.centroids(faceNo,:);
                w = w.*G.faces.areas(faceNo);
            else
                % Cell coordinates
                cellNo = rldecode((1:G.cells.num)', n, 1);
                x = cubature.transformCoords(x, cellNo, true);
                w = w.*G.cells.volumes(cellNo);
            end
        end
        
    end
    
end