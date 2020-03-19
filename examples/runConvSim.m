function output = runConvSim(G, params, varargin)
    
    opt = struct('verbose' , mrstVerbose, ...
                 'useBlock', true);
    opt = merge_options(opt, varargin{:});
    
    useBlock = opt.useBlock;
    
    force_fun = params.force_fun;
    mu_fun = params.mu_fun;
    u_fun = params.u_fun;
    eta = params.eta;
    alpha = params.alpha;
    
    doVem = false;
    [tbls, mappings] = setupStandardTables(G);
    coltbl = tbls.coltbl;
    nodefacetbl = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;
    
    Nc = G.cells.num; 
    Nf = G.faces.num; 
    Nn = G.nodes.num; 
    Nd = G.griddim; 
    
    % prepare input for analytical functions
    for idim = 1 : Nd
        cc{idim} = G.cells.centroids(:, idim);
    end
    mu = mu_fun(cc{:});
    lambda = alpha*mu;
    
    prop.mu = mu; 
    prop.lambda = lambda; 
    
    isBoundary = any(G.faces.neighbors == 0, 2); 
    bcfaces =  find(isBoundary);
    
    bcfacetbl.faces = bcfaces;
    bcfacetbl = IndexArray(bcfacetbl);
    bcnodefacetbl = crossIndexArray(bcfacetbl, nodefacetbl, {'faces'});
    bcnodefacecoltbl = crossIndexArray(bcnodefacetbl, coltbl, {}, 'optpureproduct', ...
                                       true);
    clear bcfacetbl
    
    [~, nodefacecents] = computeNodeFaceCentroids(G, tbls, eta);
    
    map = TensorMap();
    map.fromTbl = nodefacecoltbl;
    map.toTbl = bcnodefacecoltbl;
    map.mergefds = {'nodes', 'faces', 'coldim'};
    map = map.setup();
    
    bcnodefacecents = map.eval(nodefacecents);
    % Here, we assume a given structure of bcnodefacecoltbl:
    bcnodefacecents = reshape(bcnodefacecents, Nd, [])';
    bcnum = bcnodefacetbl.num;

    % Prepare input for analytical functions
    for idim = 1 : Nd
        bnfc{idim} = bcnodefacecents(:, idim);
    end
    
    % Compute boundary conditions
    for idim = 1 : Nd
        linform = zeros(bcnum, Nd);
        linform(:, idim) = 1;
        linforms{idim} = linform;
        linformvals{idim} = u_fun{idim}(bnfc{:});
    end
    
    bcfaces = bcnodefacetbl.get('faces');
    bcnodes = bcnodefacetbl.get('nodes');
    extbcnodefacetbl.faces = repmat(bcfaces, Nd, 1);
    extbcnodefacetbl.nodes = repmat(bcnodes, Nd, 1);
    extbcnodefacetbl = IndexArray(extbcnodefacetbl);
    
    bc.bcnodefacetbl = extbcnodefacetbl;
    bc.linform = vertcat(linforms{:});
    bc.linformvals = vertcat(linformvals{:});
    clear extbcnodefacetbl linforms linformvals

    % Compute body force
    force = NaN(Nc, Nd);
    for idim = 1 : Nd
        force(:, idim) = force_fun{idim}(cc{:});
    end
    % note minus sign (matter of convention)
    force = - bsxfun(@times, G.cells.volumes, force);
    % figure
    % quiver(cc{1}, cc{2}, force(:, 1), force(:, 2));
    
    % Here, we assume we know the structure of cellcoltbl;
    force = reshape(force', [], 1);
    
    loadstruct.bc = bc;
    loadstruct.force = force;
    loadstruct.extforce = zeros(tbls.nodefacecoltbl.num, 1);
    clear bc force
    
    if opt.useBlock
        assembly = blockAssembleMPSA(G, prop, loadstruct, eta, tbls, mappings, ...
                                     'blocksize', 100, 'verbose', true);
    else
        assembly = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings);
    end
    
    clear prop loadstruct
    
    B   = assembly.B  ;
    rhs = assembly.rhs;
    sol = B\rhs;

    % Displacement values at cell centers.
    cellcoltbl = tbls.cellcoltbl;
    n = cellcoltbl.num;

    u = sol(1 : n);
    u = reshape(u, Nd, [])';    

    output = struct('B'  , B  , ...
                    'rhs', rhs, ...
                    'tbls', tbls, ...
                    'u'  , u);
    
end
