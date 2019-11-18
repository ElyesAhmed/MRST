function pTop = computeTopPressure(G, p_act, flux, gamma, mu_w, k, krw)
% Get minimum pressure of the top layer using TPFA
%
% SYNOPSIS:
%   pTop = getTopPressure(G, p_act, flux, gamma, mu_w, k)
%
% PARAMETERS:
%   G       - Structure, MRST grid structure
%   p_act   - AD variable, current pressure AD-object
%   flux    - Vector, Flux corresponding to the current m-level 
%   krwAr   - Vector, Arithmetic avg of krw corresponding to current m-level
%   gamma   - Scalar, Specific gravity i.e., gamma = rho_w * g
%   mu_w    - Scalar, Dynamics viscosity
%   k       - Intrinsic permeability 
%
%  RETURNS:
%   pTop   - Scalar, minumum value of top pressure
%

%{
Copyright 2018-2019, University of Bergen.

This file is part of the fv-unsat module.

fv-unsat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fv-unsat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

nz = G.numLayers;                       % number of layers
topCellsNum = G.cells.num / nz;         % number of top cells
A  = G.faces.areas;                     % face areas
zc = G.cells.centroids(:,3);            % cell centers in z-direction 
zf = G.faces.centroids(:,3);            % face centers in z-direction
Lz = max(zf);                           % Depth of the domain
zetac = Lz - zc;                        % cell centers of elev. head
zetaf = Lz - zf;                        % face centers of elev. head
z_min = find(zf == 0);                  % idx of top faces

% Obtaining minimum value of pTop
pTop = mean( (((flux(z_min) ./ (A(z_min))) .* zc(1) .* mu_w) ./ ...
    (mean(krw(p_act.val(1:topCellsNum))) .* k)) - gamma .* ...
    (zetaf(z_min) - zetac(1:topCellsNum)) ...
    + p_act.val(1:topCellsNum));
end