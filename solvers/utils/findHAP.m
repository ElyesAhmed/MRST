function interpFace = findHAP(G, rock, bc, debug)
% Find harmonic averaging points for 2D and 3D grids. Considering both
% Dirichlet and Neumann boundary conditions

% interpFace.coords: coordinates of interpolating points
% interpFace.weights: interpolating weights
% interpFace.fraction: the fraction of cells whose centroid is
% outside the convex hull

    if nargin == 3
        debug = false;
    end

    dispif(mrstVerbose, 'findHAP\n');
    [K, ~, c] = permTensor(rock, G.griddim);
    interpFace.coords = zeros(G.faces.num, G.griddim);
    interpFace.weights = zeros(G.faces.num, 2);
    interpFace.fraction = 0;
    
    % Interior faces
    interior = all(G.faces.neighbors ~= 0, 2);
    c1 = G.faces.neighbors(interior, 1);
    c2 = G.faces.neighbors(interior, 2);
    
    n = G.faces.normals(interior, :);
    K1 = K(c1, :);
    K2 = K(c2, :);
    w1 = bsxfun(@times, K1, n(:, c));
    w2 = bsxfun(@times, K2, n(:, c));
    if G.griddim == 2
        sumfcn = @(w) [sum(w(:, 1:2), 2), sum(w(:, 3:4), 2)];
    else
        sumfcn = @(w) [sum(w(:, 1:3), 2), sum(w(:, 4:6), 2), sum(w(:, 7:9), 2)];
    end
    w1 = sumfcn(w1);
    w2 = sumfcn(w2);

    x1 = G.cells.centroids(c1, :);
    x2 = G.cells.centroids(c2, :);
    xf = G.faces.centroids(interior,:);
    xA = x1 + dot(xf-x1, n, 2) ./ dot(w1, n, 2) .* w1;
    xB = x2 + dot(xf-x2, n, 2) ./ dot(w2, n, 2) .* w2;
    w1 = vecnorm(w1, 2, 2) ./ vecnorm(xA-x1, 2, 2);
    w2 = vecnorm(w2, 2, 2) ./ vecnorm(xB-x2, 2, 2);

    interpFace.coords(interior, :) = (xA .* w1 + xB .* w2) ./ (w1 + w2);
    interpFace.weights(interior, 1) = w1;
    interpFace.weights(interior, 2) = w2;
    interpFace.weights(interior, :) = interpFace.weights(interior, :) ./ (w1 + w2);

    % Boundary faces as Hom Neumann
    bdry = ~interior;
    c1 = G.faces.neighbors(bdry, 1);
    c2 = G.faces.neighbors(bdry, 2);
    cc = max(c1, c2);
    n = G.faces.normals(bdry, :);
    K1 = K(cc, :);
    w1 = bsxfun(@times, K1, n(:, c));
    w1 = sumfcn(w1);
    x1 = G.cells.centroids(cc, :);
    xf = G.faces.centroids(bdry, :);
    xA = x1 + dot(xf-x1, n, 2) ./ dot(w1, n, 2) .* w1;
    
    interpFace.coords(bdry, :) = xA;   
    ex1 = G.faces.neighbors(:, 1) ~= 0;
    ex2 = G.faces.neighbors(:, 2) ~= 0;
    interpFace.weights(bdry & ex1, 1) = 1;
    interpFace.weights(bdry & ~ex1, 2) = 1;

        
    % interpFace2 = interpFace;


%     K = permTensor(rock, G.griddim);
%     K = reshape(K', G.griddim, G.griddim, []);
%     interpFace.coords = zeros(G.faces.num, G.griddim);
%     interpFace.weights = zeros(G.faces.num, 2);
%     interpFace.fraction = 0;

%     % Find harmonic averaging point--------------------------------------------
%     for i_face = 1:G.faces.num

%         c1 = G.faces.neighbors(i_face, 1);
%         c2 = G.faces.neighbors(i_face, 2);
%         xf = G.faces.centroids(i_face, :)';
%         if (all([c1, c2] ~= 0))
%             K1 = K(:, :, c1);
%             K2 = K(:, :, c2);
%             fn = G.faces.normals(i_face, :)';
%             w1 = K1 * fn;
%             w2 = K2 * fn;
%             x1 = G.cells.centroids(c1, :)';
%             x2 = G.cells.centroids(c2, :)';
%             xA = x1 + dot(xf-x1, fn) / dot(w1, fn) * w1;
%             xB = x2 + dot(xf-x2, fn) / dot(w2, fn) * w2;
%             w1 = norm(w1) / norm(xA-x1);
%             w2 = norm(w2) / norm(xB-x2);
%             interpFace.coords(i_face, :) = (w1 * xA + w2 * xB)' / (w1 + w2);
%             interpFace.weights(i_face, 1) = w1 / (w1 + w2);
%             interpFace.weights(i_face, 2) = w2 / (w1 + w2);
%         else
% %             ind = find(bc.face == i_face, 1);
%             %ind,bc.type{ind}
% %             if (strcmpi(bc.type{ind}, 'pressure'))
% %                 interpFace.coords(i_face, :) = xf';
% %                 interpFace.weights(i_face, (c2 == 0)+1) = bc.value{ind}(xf);
% %             else
%                 c = max(c1, c2);
%                 K1 = K(:, :, c);
%                 fn = G.faces.normals(i_face, :)';
%                 w1 = K1 * fn;
%                 x1 = G.cells.centroids(c, :)';
%                 xA = x1 + dot(xf-x1, fn) / dot(w1, fn) * w1;
%                 interpFace.coords(i_face, :) = xA';
%                 a = norm(w1) / norm(x1-xA);
%                 %gN = bc.value{ind}(xf) * G.faces.areas(i_face);
%                 gN = 0;
%                 interpFace.weights(i_face, (c1 == 0)+1) = 1;
%                 interpFace.weights(i_face, (c2 == 0)+1) = -gN / a;
%                 %interpFace.weights(i_face,(c2==0)+1)=bc.value{ind}(xf); %-gN/a;
% %             end
%         end
%     end

    % Count the number of cells whose centroid is outside the convex hull-----
    counter = zeros(G.cells.num, 1);
    for i = 1:G.cells.num
        xc = G.cells.centroids(i, :);
        theFaces = G.cells.faces(G.cells.facePos(i):G.cells.facePos(i+1)-1);
        hap = interpFace.coords(theFaces, :);
        ind = convhull(hap);
        switch G.griddim
            case 2
                xv = hap(ind, 1);
                yv = hap(ind, 2);
                counter(i) = inpolygon(xc(1), xc(2), xv, yv);
            case 3
                counter(i) = mex_inhull(xc, hap, ind, -1e-5);
        end
    end
    interpFace.fraction = 1 - sum(counter) / G.cells.num;

    if debug

        % check
        interpFace
        interpFace2
        f = interior;
        %f = bdry;
        [norm(interpFace.coords(f,:))-norm(interpFace2.coords(f,:))]
        [norm(interpFace.weights(f,:))-norm(interpFace2.weights(f,:))]
        %[norm(interpFace.coords), norm(interpFace2.coords)]
        %[norm(interpFace.weights), norm(interpFace2.weights)]

        keyboard

        interpFace.fraction
        figure
        plotGrid(G)
        hold on
        plotCellData(G, counter)
        plotHAPhull(G, interpFace, ~counter)
        keyboard
    end
end