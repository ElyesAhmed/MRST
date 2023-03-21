%% Example of poro-elastic problem
% This example sets up a poro-elastic problem which mimics a slice of the
% overburden, with an infinite horizontal well in an aquifer at the bottom
% of the domain. The poro-elastic equations are set up together with linear
% explicitely. At the end we simulate the same case by using the
% poro-elastic solver and compare with the results predicted by linear
% elasticity.
%
% More examples can be generated by modifying the parameters sent to the
% function squareTest.m, 
%

mrstModule add incomp vemmech

%% Define the grid and rock parameters
%
opt = struct( 'L'         , [10000 2000], ...
              'cartDims'  , [20 20]*1,    ...
              'grid_type' , 'square',     ...
              'disturb'   , 0.0,          ... % Parameter for disturbing grid
              'E'         , 1e9,          ... % Young's modulus
              'nu'        , 0.2);             % Poisson's ratio


G = cartGrid(opt.cartDims, opt.L);
if (opt.disturb > 0)
    G = twister(G, opt.disturb);
end
G = createAugmentedGrid(G);
G = computeGeometry(G);

Ev  = repmat(opt.E, G.cells.num, 1);
nuv = repmat(opt.nu, G.cells.num, 1);
C   = Enu2C(Ev, nuv, G);

figure()
clf, 
plotGrid(G)

%% Set up gravity as the loading term
%
gravity reset on
density = 3000; % in kg/m^3
grav    = norm(gravity());   % gravity 
load    = @(x) -(grav*density)*repmat([0, 1], size(x, 1), 1);


%% Set up the displacement boundary conditions
% 
oside = {'Left', 'Right', 'Back', 'Front'};
bc = cell(4, 1);
for i = 1 : numel(oside);
    bc{i} = pside([], G, oside{i}, 0);
    bc{i} = rmfield(bc{i}, 'type');
    bc{i} = rmfield(bc{i}, 'sat');
end
% Find the nodes for the different sides and set the boundaray conditions for
% elastisity.
for i = 1 : 4
    inodes = mcolon(G.faces.nodePos(bc{i}.face), G.faces.nodePos(bc{i}.face + 1) - 1);
    nodes = unique(G.faces.nodes(inodes));
    disp_bc = struct('nodes'   , nodes,      ...
                     'uu'      , 0,          ...
                     'faces'   , bc{i}.face, ...
                     'uu_face' , 0,          ...
                     'mask'    , true(numel(nodes), G.griddim));
    bc{i}.el_bc = struct('disp_bc', disp_bc, 'force_bc', []);
end
bcdisp = @(x) x*0.0; % Boundary displacement function set to zero.
% Setup diriclet boundary conditions at selected sides
%
% On the left and right-hand sides, zero displacement is imposed in the x
% direction, and free in the y direction. This is done by using masks.
bc_el_sides{1} = bc{1}; 
bc_el_sides{1}.el_bc.disp_bc.mask(:, 2) = false;
bc_el_sides{2} = bc{2};
bc_el_sides{2}.el_bc.disp_bc.mask(:, 2) = false;
% Zero displacement is imposed in the y direction at the bottom, free
% displacement at the top.
bc_el_sides{3} = bc{3};
bc_el_sides{3}.el_bc.disp_bc.mask(:, 1) = false;
bc_el_sides{4} = [];

% Collect the displacement boundary conditions
nodes = [];
faces = [];
mask = [];
for i = 1 : numel(bc)
    if(~isempty(bc_el_sides{i}))
        nodes = [nodes; bc_el_sides{i}.el_bc.disp_bc.nodes]; %#ok
        faces = [faces; bc_el_sides{i}.el_bc.disp_bc.faces]; %#ok
        mask  = [mask; bc_el_sides{i}.el_bc.disp_bc.mask]; %#ok
    end
end
disp_node = bcdisp(G.nodes.coords(nodes, :));
disp_faces = bcdisp(G.faces.centroids(faces, :));
disp_bc = struct('nodes', nodes, 'uu', disp_node, 'faces', faces, 'uu_face', disp_faces, 'mask', mask); 

%% Set up the force boundary conditions

% A force is applied on the top surface. 
sigma = opt.L(2)/10;
force = 100*barsa;
face_force = @(x) force*exp(-(((x(:, 1) - opt.L(1)/2))./sigma).^2) + 10*barsa;
faces = bc{4}.face;
% Set up the force boundary structure 
% Note that the unit for the force is in units Pa/m^3
force_bc = struct('faces', faces, 'force', bsxfun(@times, face_force(G.faces.centroids(faces, :)), [0 -1]));    


% Final structure for the  boundary conditions
el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);


%% Assemble and solve the system
%
[uu, extra] = VEM_linElast(G, C, el_bc, load);

% We retrieve the discrete operators that have been assembled
As        = extra.disc.A;     % Global stiffness matrix
div       = extra.disc.div;   % Divergence operator
gradP     = extra.disc.gradP; % Gradient operator (dual of div)
isdirdofs = extra.disc.isdirdofs; % Degrees of freedom where in fact
                                  % Dirichlet boundary conditions are imposed.
rhs_s     = extra.disc.rhs;   % Right-hand side
Vdir      = extra.disc.V_dir; %
ind_s     = [1 : size(As, 1)]';%#ok


%% Set up the parameters for the flow part
%
perm = 1*darcy*ones(G.cartDims);
perm(:, floor(G.cartDims(2)/5) : end) = 1*milli*darcy;
rock = struct('perm', reshape(perm, [], 1),     ...
              'poro', 0.1*ones(G.cells.num, 1), ...
              'alpha', ones(G.cells.num, 1)); % Rock parameters
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1000); % Fluid parameters
rock.cr = 1e-4/barsa; % Fluid compressibility
T = computeTrans(G, rock); % Compute the transmissibilities
pv = poreVolume(G, rock); % Compute the pore volume
pressure = 100*barsa*ones(G.cells.num, 1); % Initial pressure
% Set up the initial state for the flow. 
state = struct('pressure', pressure, 's', ones(G.cells.num, 1), 'flux', zeros(G.faces.num, 1));
% Set up a well structure for the fluid injection and also boundary
% conditions for the flow.
mcoord = [5000 200];
[dd, wc] = min(sum(bsxfun(@minus, G.cells.centroids, mcoord).^2, 2));
W = addWell([], G, rock, wc, 'type', 'bhp', 'val', 3000*barsa);
bc_f = pside([], G, 'Front', 10*barsa);
fbc = addFluidContribMechVEM(G, bc_f, rock, isdirdofs);

%% Set up matrix for mechanics and flow
% Use a TPFA solver for weakly compressible flow. From the output we
% retrieve the operators for the fluid.
dt = day; % time step.
state = lincompTPFA(dt, state, G, T, pv, fluid, rock, 'MatrixOutput', true, 'wells', W, 'bc', bc_f);

Af    = state.A;   % Global matrix for the flow solver
orhsf = state.rhs; % Right-hand side for the flow equations (includes well)
ct    = state.ct;  % Compressibility coefficient
ind_f = [ind_s(end) + 1 : ind_s(end) + G.cells.num]';%#ok

% x denotes the global unknown (displacement + pressure)
x = zeros(ind_f(end), 1);
x(ind_f) = pressure;
p = pressure;

uu0 = uu;            % Solution from linear elasticity
u_tmp = reshape(uu', [], 1);
x(1 : ind_s(end)) = u_tmp(~isdirdofs);
u = zeros(numel(isdirdofs), 1);
rhsf = zeros(size(orhsf));
plotops = {'EdgeColor', 'none'};
fac = rock.poro(1); 

% Assemble the matrix for the global system, mechanics + flow
zeromat = sparse(size(Af, 1) - G.cells.num, size(div, 2));
SS = [As, [fac*(-gradP), zeromat']; [fac*div; zeromat], ct + dt*Af];

%% Time loop
% Set up bigger figures
df = get(0, 'DefaultFigurePosition');
figure(1); set(1, 'Position', df.*[0.8, 1, 3, 1]); clf
figure(2); set(2, 'Position', df.*[1, 1, 3, 1]); clf

t = 0;
end_time = 10*dt;
while t < end_time
    t = t + dt;
    % We update the right-hand side for the fluid (divergence term comes
    % from mechanics).
    rhsf(1 : G.cells.num) = orhsf(1 : G.cells.num)*dt + ct(1 : G.cells.num, 1 : G.cells.num)*p + fac*div*x(ind_s);
    rhsf(G.cells.num + 1 : end) = orhsf(G.cells.num + 1 : end);
    rhs = [rhs_s - fbc; rhsf];

    % Solve the system
    x = SS\rhs;

    % Retrieve pressure and displacement
    p = x(ind_f); % retrieve pressure
    u(isdirdofs) = Vdir(isdirdofs); 
    u(~isdirdofs) = x(ind_s); % retrieve displacement
    uu = reshape(u, G.griddim, [])';

    % Plot absolute results from poromechanics
    figure(1)
    subplot(2, 2, 1), cla
    plotNodeData(G, uu(:, 1), plotops{ : }); colorbar;
    title('Displacement in x direction')
    subplot(2, 2, 2), cla
    plotNodeData(G, uu(:, 2), plotops{ : }); colorbar;
    title('Displacement in y direction')
    subplot(2, 2, 3), cla
    plotCellDataDeformed(G, p/barsa, uu); colorbar;
    title('Pressure')
    subplot(2, 2, 4), cla
    ovdiv = extra.disc.ovol_div;
    mdiv = ovdiv*reshape(uu', [], 1)./G.cells.volumes;
    plotCellDataDeformed(G, mdiv, uu); colorbar();
    title('Divergence')
    
    % Plot difference from results predicted by linear elasticity
    figure(2)
    uur = uu - uu0;
    subplot(2, 2, 1), cla
    plotNodeData(G, uur(:, 1), plotops{ : }); colorbar;
    title('Relative displacement in x direction')
    subplot(2, 2, 2), cla
    plotNodeData(G, uur(:, 2), plotops{ : }); colorbar;
    title('Relative displacement in y direction')
    subplot(2, 2, 3), cla
    plotCellDataDeformed(G, p/barsa, uur); colorbar;
    title('Pressure (grid deformed using relative disp.)')
    axis tight
    subplot(2, 2, 4), cla
    ovdiv = extra.disc.ovol_div;
    mdiv = ovdiv*reshape(uur', [], 1)./G.cells.volumes;
    plotCellDataDeformed(G, mdiv, uu); colorbar();   
    title('Divergence of relative displacement')
    axis tight
    %pause(0.01);
    
end

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
