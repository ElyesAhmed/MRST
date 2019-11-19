function p0 = pointsSingleWellNode(pW, ly, ny, na, ii)
% Generate the 2D points corresponding to single well node ii
% The parameter represeation can be found in 'prepareWellRegionNodes2D'
    
    % Compute the angle
    y0 = linspace(ly/2, 0, ny/2+1)';
    y0 = y0(1:end-1);
    
    fTheta = @(x,y)2*pi*double(sign(atan2(y,x))<0) + atan2(y,x);
    pW = pW(:, [1,2]);
    p2 = pW(ii, :);
    if ii == 1
        p3  = pW(ii+1, :);
        ang = fTheta(p3(1)-p2(1), p3(2)-p2(2));
        g1 = 0.5*pi;
        g2 = 1.5*pi;
    elseif ii == size(pW,1)
        p1  = pW(ii-1, :);
        ang = fTheta(p2(1)-p1(1), p2(2)-p1(2));
        g1 =  0.5*pi;
        g2 = -0.5*pi;
    else
%         p3   = pW(ii+1, :);
%         ang = fTheta(p3(1)-p2(1), p3(2)-p2(2));
        p1  = pW(ii-1, :);
        p3  = pW(ii+1, :);
        p12 = (p1+p2)/2;
        p23 = (p2+p3)/2;
        ang = fTheta(p23(1)-p12(1), p23(2)-p12(2));
    end

    % Rotating mapping
    M = [cos(ang), sin(ang); ...
        -sin(ang), cos(ang)];

    % Cartesian points
    y  = [y0; 0; -y0(end:-1:1)];
    xy = [zeros(size(y)), y];
    xy = xy * M;
    xy = bsxfun(@plus, xy, p2);
    
    % Radial points
    if ii == 1 || ii == size(pW,1)
        g = linspace(g1, g2, na+1);
        g = g(2:end-1);
        r = y0;
        xr  = r * cos(g); 
        yr  = r * sin(g);
        xyr = [xr(:), yr(:)];
        xyr = xyr * M;
        xyr = bsxfun(@plus, xyr, p2);
    else
        xyr = [];
    end

    p0.cart = xy;
    p0.rad  = xyr;
end