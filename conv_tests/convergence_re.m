function [errorPsi, errorFlux] = convergence_re(cells, timeLevels)
% Converge of Richards' equation (Numerical Example 1 from the Chapter)
%
% SYNOPSIS:
%   function [error_psi, error_flux] = convergence_re(cells, timeLevels)
%
% PARAMETERS:
%   cells       - Scalar, number of cells = nx = ny.
%   timeLevels  - Scalar, number of discrete time levels.
%
%  RETURNS:
%   errorPsi   - Scalar, L2-like discrete error for the pressure.
%   errorFlux  - Scalar, L2-like discrete error for the flux.
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

% Importing modules
mrstModule add fvbiot

% Creating the grid
nx = cells;                         % Cells in x-direction
ny = cells;                         % Cells in y-direction
Lx = 1;                             % Lenght in x-direction
Ly = 1;                             % Length in y-direction
G = cartGrid([nx, ny], [Lx, Ly]);   % Create Cartesian grid
G = computeGeometry(G);             % compute geometry

% Extracting grid information
V = G.cells.volumes;                % Cell volumes
A = G.faces.areas;                  % Face areas
xc = G.cells.centroids(:, 1);       % cell centers in x-direction
yc = G.cells.centroids(:, 2);       % cell centers in y-direction
xf = G.faces.centroids(:, 1);       % face centers in x-direction
yf = G.faces.centroids(:, 2);       % face centers in y-direction
nu_x = G.faces.normals(:, 1);       % face normals in x-direction
nu_y = G.faces.normals(:, 2);       % face normals in y-direction

% Physical parameters

% Rock
k = 1;  % intrinsic permeability
rock.perm = repmat([k, k], [G.cells.num, 1]); % creating perm structure

% Water retention curves (van Genuchten - Mualem)
alpha = 0.04;       % Equation parameter
nVan = 2;           % Equation parameter
mVan = 1-(1/nVan);  % Equation parameter
theta_s = 0.4;      % Water content at saturation conditions
theta_r = 0.1;      % Residual water content 

%%Water retention curves

[theta, krw, C_theta] = vanGenuchtenMualemTheta(alpha, theta_s, ...
                                                theta_r, nVan, mVan);

% Boundary and Initial Conditions

% Boundary indices
x_min = find(xf == 0);                           % west faces
x_max = find(xf > 0.9999*Lx & xf < 1.0001*Lx );  % east faces

y_min = find(yf == 0);                           % south faces 
y_max = find(yf > 0.9999*Ly & yf < 1.0001*Ly );  % north faces 

% Boundary cond. structure
boundPsi = -1;
bc = addBC([], x_min, 'pressure', boundPsi);  
bc = addBC(bc, x_max, 'pressure', boundPsi);  
bc = addBC(bc, y_min, 'pressure', boundPsi);  
bc = addBC(bc, y_max, 'pressure', boundPsi);  

% Boundary cond. values 
bc_val = zeros(G.faces.num, 1);   % initializing
bc_val(x_min) = boundPsi;  % west faces
bc_val(x_max) = boundPsi;  % east faces
bc_val(y_min) = boundPsi;  % south faces
bc_val(y_max) = boundPsi;  % north faces

% Initial condition
psi_init = boundPsi * ones(G.cells.num, 1);       

% Calling MPFA routine and creating operators

mpfa_discr = mpfa(G, rock, [], 'bc', bc, 'invertBlocks', 'matlab');
     
F = @(x) mpfa_discr.F * x;               % flux discretization
boundF = @(x) mpfa_discr.boundFlux * x;  % boundary discretization
divF = @(x) mpfa_discr.div * x;          % divergence of the flux

% Creating AD variable
psi_ad = initVariablesADI(psi_init);

% Time parameters
iniTime = 0;  % intial simulation time
simTime = 1;  % final simulation time 
times = linspace(iniTime, simTime, timeLevels + 1); % evaluation times
tau = diff(times); % time steps

% Analytical solution + Source term

% Note that here we make use of the Matlab Symbolic Toolbox.
% This is not ideal, but is more practical than hard-coding a very
% long function. It also allows the user to make changes in the analytical
% solution to determine convergence rates for different scenarios.

syms x y t

% Analytical solution
psi = -t * x * (1 - x) * sin(pi*x) * y * (1 - y) * cos(pi*y) + boundPsi;

% Parameters
wc = (theta_s - theta_r) / (1 + (alpha * abs(psi)) ^ nVan) ...
     ^ mVan + theta_s;
 
rp = (1 - (alpha * abs(psi)) ^ (nVan-1) ...
     * (1 + (alpha * abs(psi)) ^ nVan) ^ (-mVan)) ^ 2 ...
     / (1 + (alpha * abs(psi)) ^ nVan) ^ (mVan / 2);

% Obtaining analytical expressions
grad_psi = [diff(psi, x); diff(psi, y)];
vel = -rp * grad_psi;
div_vel = diff(vel(1), x) + diff(vel(2), y);
partial_theta_partial_t = diff(wc, t);
source = partial_theta_partial_t + div_vel;

% Converting into inline functions
psi_ex = matlabFunction(psi);
f_ex = matlabFunction(source);
q_ex = matlabFunction(vel);

% Discrete equations
                          
% Arithmetic average
krwAr = @(psi_m) arithmeticAverageMPFA(G, bc, krw, psi_m);      

% Darcy Flux
q = @(psi, psi_m)  krwAr(psi_m) .* (F(psi) + boundF(bc_val));

% Mass Conservation                         
psiEq = @(psi, psi_n, psi_m, tau, source) (V ./ tau) .* ( ... 
                theta(psi_m) + C_theta(psi_m) .* (psi - psi_m) - theta(psi_n) ...
            ) + divF(q(psi, psi_m)) - V .* source;
                            
% Time loop

tt = 1;         % time counter
tol = 1E-8;     % tolerance
maxIter = 10;   % maximum number of iterations

while times(tt) < simTime
    
    psi_n = psi_ad.val; % current time level (n-index)    
    tt = tt + 1;        % increasing time counter
    source = f_ex(times(tt), xc, yc);  % Obtaining source term

    % Calling newton solver
    [psi_ad, psi_m, ~] = solverRE(psi_ad, psi_n, psiEq, ...
                                        tau(tt-1), source, times(tt), ...
                                        tol, maxIter);                                                                                 
end

% Collecting results and computing error

% True solution
psi_true = psi_ex(simTime, xc, yc);           % Exact pressure head
q_true = q_ex(simTime, xf, yf);               % Exact velocity field
Q_true = q_true(1:G.faces.num) .* nu_x ...    % Exact (normal) fluxes
         + q_true(G.faces.num+1:end) .* nu_y; %

% Numerical solution
psi_num = psi_ad.val;                         % Numerical pressure head
Q_num = q(psi_ad.val, psi_m);                 % Numerical fluxes

% Computing errors
errorPsi = sqrt(sum(V .* (psi_true - psi_num).^2)) ... 
            ./ sqrt(sum(V .* psi_true.^2));
errorFlux = sqrt(sum(A .* (Q_true - Q_num).^2)) ... 
            ./ sqrt(sum(A .* Q_true.^2));
