function G = computeCellDimensions(G)

    % Get extra grid topology
    G = createAugmentedGrid(G);

    % Get node coordinates
    xn  = G.nodes.coords(G.cells.nodes,:);
    ncn = diff(G.cells.nodePos);
    
    % Get minimum and maximum cell coordinates
    [G.cells.xMin, G.cells.xMax] = getMinMax(xn, ncn);

    dx = max(G.cells.xMax - G.cells.centroids, G.cells.centroids - G.cells.xMin)*2;
    % Compute cell bounding box dimensions
    G.cells.dx = dx;%G.cells.xMax - G.cells.xMin;
    G.cells.basisCenters = G.cells.xMin + G.cells.dx/2;
    
    % Get face coordinates
    xn  = G.nodes.coords(G.faces.nodes,:);
    nfn = diff(G.faces.nodePos);
    xn = xn - rldecode(G.faces.centroids, nfn, 1);
    if G.griddim == 3
        % Create local coordinate systems on each face
        G.faces.coordSys = faceCoordSys(G);
        
        % Map to face coordinate system
        node2face = rldecode((1:G.faces.num)', nfn, 1);
        xnf       = zeros(sum(nfn), 2); 
        for dNo = 1:2
            vec = G.faces.coordSys{dNo}(node2face,:);
            xnf(:, dNo) = sum(xn.*vec,2);
        end
        xn = xnf;
        
    end
        
    % Get minimum and maximum cell coordinates
    [G.faces.xMin, G.faces.xMax] = getMinMax(xn, nfn);
    G.faces.dx = G.faces.xMax - G.faces.xMin;
    
    xbr = G.faces.xMin + G.faces.dx/2;
    xb = 0;
    if G.griddim == 3
        for dNo = 1:2
            vec = G.faces.coordSys{dNo};
            xb = xb + xbr(:,dNo).*vec;
        end
    end
    xb = xb + G.faces.centroids;
    G.faces.basisCenters = xb;
    
    G.faces.phys2ref = @(x,faces) phys2ref(x,G,faces);
    G.faces.ref2phys = @(x,faces) ref2phys(x,G,faces);
        
    if 0
        for f = 1:G.faces.num
            clf;
            plotFaces(G, f);
            hold on
            x = G.faces.basisCenters(f,:);
            plot3(x(:,1), x(:,2), x(:,3), '.');
            axis equal tight
        end
    end
        
end

function coordSys = faceCoordSys(G)
    
    % Get nodes of last edge of each face
    nfe   = diff(G.faces.edgePos);
    edges = G.faces.edges(cumsum(nfe));
    nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), ...
                                 G.edges.nodePos(edges+1)-1));

    % First coordinate axis points along last edge
    vec1 = G.nodes.coords(nodes(2:2:end),:) ...
         - G.nodes.coords(nodes(1:2:end-1),:);
    vec1 = vec1./sqrt(sum(vec1.^2,2));
    
    % Second vector constructed as cross product of vec1 and face normal
    n    = G.faces.normals./G.faces.areas;
    vec2 = cross(vec1, n, 2);
    
    coordSys = {vec1, vec2};

end

function xr = phys2ref(x, G, faces)

    x  = x - G.faces.basisCenters(faces,:);
    xr = zeros(numel(faces),2);
    for dNo = 1:2
        xr(:,dNo) = sum(x.*G.faces.coordSys{dNo}(faces,:),2);
    end
    xr = xr./(G.faces.dx(faces,:)/2);

end

function x = ref2phys(xr, G, faces)

    xr = xr.*G.faces.dx(faces,:)/2;
    x  = 0;
    for dNo = 1:2
        x = x + G.faces.coordSys{dNo}(faces,:).*xr(:,dNo);
    end
    x = x + G.faces.basisCenters(faces,:);

end
