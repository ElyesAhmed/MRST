function grdecl = convertAtlasTo3D(meta_thick, meta_top, data_thick, data_top, nz)
%Create GRDECL struct from CO2 storage atlas thickness/top data
%
% SYNOPSIS:
%   grdecl = convertAtlasTo3D(m_thick, m_top, d_thick, d_top, 3)
%
% DESCRIPTION:
%   Given two datasets with possibily non-matching nodes, interpolate and
%   combine to form a GRDECL struct suitable for 3D simulations after a
%   call to processGRDECL or for further manipulation.
%
% REQUIRED PARAMETERS:
%   meta_thick - Metainformation for the thickness data as produced from
%                readAAIGrid
%
%   meta_top   - Metainformation for top data.
%
%   data_thick - Data for the thickness data as produced from
%                readAAIGrid
%
%   top_thick  - Data for the top data as produced from
%                readAAIGrid
%
%
% RETURNS:
%   grdecl - GRDECL struct suitable for processGRDECL
%
%
% NOTES:
%   It is likely easier to use getAtlasGrid instead of calling this
%   routine directly.
%
% SEE ALSO:
%   getAlasGrid

%{
#COPYRIGHT#
%}

    ndims = [meta_top.dims, nz + 1];
    dims = ndims - 1;
    h = meta_top.cellsize;

    
    grdecl.cartDims = reshape(dims, 1, []);
    
    xl = meta_top.xllcorner;
    yl = meta_top.yllcorner;
    
    % Create grids
    [X, Y, Z]  = ndgrid(linspace(xl, dims(1)*h + xl, ndims(1)), ...
                        linspace(yl, dims(2)*h + yl, ndims(2)), ...
                        linspace(0,  1,              ndims(3)));
    

    F_top   = interpolateData(meta_top,   data_top);
    F_thick = interpolateData(meta_thick, data_thick);
    %if(false)
    %   x = squeeze(X(:,:,1));
    %   y = squeeze(Y(:,:,1));  
    %   thick = F_thick(x, y)/dims(3);
    %   top   = F_top(x,y);
    %else
    x = reshape(X(:,:,1), [], 1);
    y = reshape(Y(:,:,1), [], 1);
    thick = reshape(F_thick(x, y), ndims(1), ndims(2))/dims(3);
    top   = reshape(F_top(x, y)  , ndims(1), ndims(2));
    %end
    
   
    
    %
    
    thick(thick<=0) = NaN;
    %
    
    % Uniformly partition the thickness across the layers
    for i = 1:ndims(3)
        Z(:,:,i) = top + thick*(i-1);
    end
    Z = sort(Z, 3);
    
    
    % Make pillars
    n = prod(ndims(1:2));
    lines = zeros([n, 6]);
    lines(:, [1, 4]) = reshape(X(:,:,[1, end]), [n, 2]);
    lines(:, [2, 5]) = reshape(Y(:,:,[1, end]), [n, 2]);
    lines(:, [3, 6]) = reshape(Z(:,:,[1, end]), [n, 2]);
    grdecl.COORD = reshape(lines.', [], 1);

    % Assign z-coordinates
    % ind(d) == [1, 2, 2, 3, 3, ..., dims(d), dims(d), dims(d)+1]
    ind = @(d) 1 + fix((1 : 2*dims(d)) ./ 2);
    z   = Z(ind(1), ind(2), ind(3));

    grdecl.ZCORN = z(:);

    % Assign active cells by removing cells containing NaN points
    
    z = squeeze(isnan(Z(:,:,1)));
    
    tmp = z;
    for i = -1:1
        for j = -1:1
            % mask away any cells touching nan-points
            tmp = tmp | circshift(z, [i,j]);
        end
    end
    
    tmp = ~(tmp(1:end-1, 1:end-1));
    for i = 2:dims(3)
        tmp(:,:,i) = tmp(:,:,1);
    end
    grdecl.ACTNUM = reshape(int32(tmp), [], 1);
    
    grdecl.ZCORN(isnan(grdecl.ZCORN)) = inf;
    grdecl.COORD(isnan(grdecl.COORD)) = inf;
end

function F = interpolateData(meta, data)
    dims = meta.dims;
    % We have a cell centered grid, so subtract by one
    gdims = dims - 1;
    h = meta.cellsize;
    xl = meta.xllcorner;
    yl = meta.yllcorner;
    [X, Y]  = ndgrid(linspace(xl, gdims(1)*h + xl, dims(1)), ...
                     linspace(yl, gdims(2)*h + yl, dims(2)));    
    F =@(x,y) interp2(X',Y',data', x, y);              
    %F = TriScatteredInterp(X(:), Y(:), data(:));
end
