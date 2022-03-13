function [Res, Jac] = twophaseJacobianWithVE_coupled(G, state, rock, fluid, varargin)
%Residual and Jacobian of single point upwind solver for two-phase flow.
%
% SYNOPSIS:
%   resSol = twophaseJacobian(G, state, rock, fluid)
%   resSol = twophaseJacobian(G, state, rock, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophaseJacobian returns function handles for the residual
%   and its Jacobian matrix for the implicit upwind-mobility weighted
%   dicretization of
%              __
%        s_t + \/· [f(s)(v·n + mo(rho_w - rho_o)n·Kg)] = f(s)q
%%
%   where v·n is the sum of the phase Dary fluxes, f is the fractional
%   flow function,
%
%                  mw(s)
%        f(s) = -------------
%               mw(s) + mo(s)
%
%   mi = kr_i/mu_i is the phase mobiliy of phase i, mu_i and rho_i are the
%   phase viscosity and density, respectivelym, g the (vector) acceleration
%   of gravity and K the permeability.  The source term f(s)q is a
%   volumetric rate of water.
%
%   Using a first-order upwind discretisation in space and a backward Euler
%   discretisation in time,  the residual of the nonlinear system of
%   equations that must be solved to move the solution state.s from time=0
%   to time=tf, are obtained by calling F(state, s0, dt) which yields
%   Likewise, the Jacobian matrix is obtained using the function Jac.
%
% REQUIRED PARAMETERS:
%   resSol  - Reservoir solution structure containing valid
%             saturation resSol.s with one value for each cell in
%             the grid.
%
%   G       - Grid data structure discretising the reservoir model.
%
%   rock    - Struct with fields perm and poro.
%
%   fluid   - Data structure describing the fluids in the problem.
%
% OPTIONAL PARAMETERS:
%   - verbose 
%
%   - wells
%
%   - src
%
%   - bc
%
%   - vert_avrg : if true use vertical average formulation of gravity and
%                  capillary forces, need a suitable fluid object
%                  default false
%
%   - vert_method : method used for vertical average on 3d grids
%                   valid options are 
%                  'topface': use top surface for gravity gradient
%                  'cells' : use cellcentroid to
%                  'pp_cells': use cellcentroid just a bit different
%
%
% RETURNS:
%   F       - Residual
%
%   Jac     - Jacobian matrix (with respect to s) of residual.
%
% EXAMPLE:
%
%
%
% SEE ALSO:
%   `implicitTransport`, `explicitTransport`.

% TODO:
%   - implement gravity effects for pressure boundary and wells
%   - multipliers for gravity flux.  Handle this in 'getFlux'

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

% $Date: 2012-10-05 10:15:39 +0200 (Fri, 05 Oct 2012) $
% $Revision: 10006 $

% NOTE we solve
%              __
%        s_t + \/· [D·f(s) + D·f(s)·mo(s)] = f(s) Q
%
% for constant vector and scalar fields fields D, G and Q.


   opt = struct('verbose', false,... 
                 'wells', [],...
                 'src', [],...
                 'bc', [],...
                 'vert_avrg',false,...
                 'vert_method','topface', ...
                 'region3D', false(G.cells.num,1));
   opt = merge_options(opt, varargin{:});

   assert (size(state.s,2)<3 || all(state.s(:,3)==0));


   % All source terms, i.e., boundary conditions, source terms.
   % Wells are handeled explicitly further down.
   q = computeTransportSourceTerm(state, G, [], opt.src, opt.bc);
   q = q(:);

   % Extract (constant) fluid densities.
   [rho, rho]     = fluid.properties(state);

   % Compute (constant) parts of flux, i.e., the total Darcy velocity v·n
   % and the gravitational potential (rho_w - rho_o)n·Kg across each
   % internal face in the grid.
   if ~isfield(fluid, {'pc'})
      [dflux, gflux] = getFlux(G, rock, rho, state.flux);
      pcflux = @(s)( zeros(size(gflux)));
      pcJac=0;
     
      % Bind flux terms dflux and gflux.  The call [iw,io]=findUpwindCells(..)
      % computes cell index of upwind cell for the water phase (iw) and oil
      % phase (io).
      cdata           = computeConstData(G, dflux, gflux, opt);
      findUpwindCells = @(mob, resSol) findPhaseUpwindCells(cdata, mob);      
   else
      [dflux, gflux, pcflux, pcJac] = getFluxCap(G,state, rock, rho,...
                                                  fluid, opt.vert_avrg,opt.vert_method, opt.region3D);
       % take into account the change in capillary forces in each step..
      findUpwindCells = @(mob, resSol) findPhaseUpwindCells(computeConstData(G, ...
                                    dflux, gflux+pcflux(resSol), opt), mob);
   end

    pv = poreVolume(G, rock);


   % Return function handles
   Res = @(xr, xr0, dt) Residual(xr, xr0, dt, fluid, dflux, gflux, pcflux,...
                                findUpwindCells, G, pv, q, opt.wells);

   Jac = @(xr, xr0, dt)     Jacobian(xr, xr0, dt,     fluid, dflux, gflux, pcflux, pcJac,...
                                findUpwindCells, G, pv, q, opt.wells);
end


%--------------------------------------------------------------------------
function J = Jacobian (resSol,resSol_0, dt, fluid, dflux, gflux, pcflux, pcJac, findIx, G, pv, q, W)
%  At an interface, we compute the Jacobain of F as follows,
%

%           dt    mw
%     F =   -- · ------- (dflux + mo*(pcflux + gflux))
%           pv   mw + mo 
%

%    dF      dt      mo       1
%    --- =   -- · -------  -----  (dflux + mo*(pcflux + gflux))
%    dmw     pv   mw + mo  mw + mo
%
%            dt  fo
%        =   --  --  (dflux + mo*(pcflux + gflux))
%            pv  mt
%
%
%    dF      dt  fw
%    --- = - --  --  (dflux - mw*(pcflux + gflux))
%    dmo     pv  mt
%
%
%    dF          dt  fw
%    ------- =   --  --  mo
%    dpcflux     pv  
%
%
%    dF    dF  dmw    dF  dmo   dF      dpcflux
%    --  = --- ---  + --- --- + ------- -------
%    dS    dmw dS     dmo dS    dpcflux   dS
%
%  If the mobilities mw and mo are evaluated in grid cells iw and io,
%  respectively, then the respective (matrix) contributions of each
%  term is added to elements (i, iw) and (i, io) of the Jacobian matrix.
    if isfield(resSol, 'extSat'),
         % Save extremal values of saturation for hysteresis purposes
         % First column contains minima, second maxima
           resSol.extSat(:,1) = min(resSol.s(:,1), resSol_0.extSat(:,1));
           resSol.extSat(:,2) = max(resSol.s(:,1), resSol_0.extSat(:,2));
   end
   mu        = fluid.properties(resSol);     % viscosity
   sat       = fluid.saturation(resSol);     % cell saturation
   [kr, dkr] = fluid.relperm(sat, resSol);
   m         = bsxfun(@rdivide, kr,  mu);    % mobility in each cell
   dm        = bsxfun(@rdivide, dkr, mu);    % derivative of mobility
   clear mu sat kr dkr

   internal_faces   = all(G.faces.neighbors~=0, 2);
   ic1  = double(G.faces.neighbors(internal_faces,1));
   ic2  = double(G.faces.neighbors(internal_faces,2));
   
   % For each face k, iw(k) and io(k) is the cell index in which we
   % must evaluate the water and oil mobilities, respectively.
   % NB: always use 3D cell for upwind on 2D/3D boundary
   [iw, io] = findIx(m, resSol);
   
   
   
%   if(isfield(resSol,'p_ph'))      
%      [iww,ioo]=upwindPh(G,resSol,rho);
%   end
   % % temporary removed:
%    if ~any(gflux),
%       assert(all(iw == io));
%    end

   m_face   = [m( iw, 1) m(io,2)];
   dm_face  = [dm(iw,1) dm(io,2)];
   [m_face dm_face] = fixCoupledBoundaryMob(G, m_face, resSol, fluid, iw, io, dm_face);
   f_face   = bsxfun(@rdivide, m_face, sum(m_face,2));

   dmw_fo = f_face(:,2).*dm_face(:,1);
   dmo_fw = f_face(:,1).*dm_face(:,2);

   dvw_ds1 =  (dflux + (pcflux(resSol) + gflux).*m_face(:,2)) .* dmw_fo ./ sum(m_face,2);
   dvw_ds2 = -(dflux - (pcflux(resSol) + gflux).*m_face(:,1)) .* dmo_fw ./ sum(m_face,2);

   
   dFdSiw1  = dt./pv(ic1).*dvw_ds1;
   dFdSiw2  = dt./pv(ic2).*dvw_ds1;
   dFdSio1  = dt./pv(ic1).*dvw_ds2;
   dFdSio2  = dt./pv(ic2).*dvw_ds2;

      
   clear dmw_fo dmo_fw dvw_ds1 dvw_ds2 internal_faces

   % In cell: differentiate term s-s0 + dt/pv(max(q,0) + min(q,0).*f) with
   % respect to cell saturation.
   f   = bsxfun(@rdivide, m, sum(m,2));
   df  = (f(:,2).*dm(:,1) - f(:,1).*dm(:,2))./sum(m,2);
   
      
   % derivative of Sn+1 and contribution from in perforations of well
   d   = 1-dt./pv.*(min(q, 0).*df);

   id  = (1:G.cells.num)';

   % Assemble contribututions 
   %           diagonal  water mobility      oil mobility
   J = sparse([id;        ic1;      ic2;        ic1;      ic2], ... %row
              [id;         iw;       iw;         io;       io], ... %column
              [d;     dFdSiw1; -dFdSiw2;    dFdSio1; -dFdSio2], ... %value
              G.cells.num, G.cells.num);
   
   
   % Add contributions from well for cases where there is inflow (and
   % possibly also outflow) in the well)
   % Calculate part of jacobian assosiated with wells
   iout=[];iinn=[];matel=[];
   for kk=1:numel(W)
       inperf  =  resSol.wellSol(kk).flux<0;
       outperf =  resSol.wellSol(kk).flux>0;
       if(sum(inperf)>0) 
          %S_avg = -resSol.wellSol(kk).flux(inperf).*f(W(kk).cells(inperf));%+W(kk).comp(1);
          influx = -sum(resSol.wellSol(kk).flux(inperf)); % sum influx (positive) 
          wellhead = sum(resSol.wellSol(kk).flux);        % net injection or prodcution  
          
          if(wellhead>0) % well is net injector (more outflux than influx)
             normal=(influx+wellhead); % outflux
             %S_avg = (sum(S_avg)+W(kk).val*W(kk).compi(1))./(normal);
          else % well is net producer
             normal=influx; % influx
             %S_avg = (sum(S_avg))./(normal);
          end
                  
          flux_in  = resSol.wellSol(kk).flux(inperf);
          flux_out = resSol.wellSol(kk).flux(outperf);    
          
          num_in  = sum(inperf);
          num_out = sum(outperf); 
          
          % compute source term
          q(W(kk).cells(inperf)) = q(W(kk).cells(inperf)) + flux_in;
         
          % If the well also has outflow cells, we must add 
          % Jacobi contributions to out cells from the in cells
          if(sum(outperf)>0) 
             % row
             iout = [iout;rldecode(W(kk).cells(outperf),repmat(num_in, num_out,1))];
             % column
             iinn = [iinn; repmat(W(kk).cells(inperf),num_out,1)];
             
             % each in-perforation gives contribution to all out cells. The
             % contribution for one out cell is weigthed by the relative
             % outflow from the well into that cell given by: 
             val_tmp = bsxfun(@times, repmat(flux_out,1, num_in ), ...
                           ((-flux_in.*df(W(kk).cells(inperf)))/normal)');                        
             % value
             matel  =[matel; reshape(val_tmp', [], 1) ];
          end        
       end
       %q(W(kk).cells(outperf)) = q(W(kk).cells(outperf))
   end
   Jwell = sparse(iout,iinn,-matel,G.cells.num,G.cells.num);
   
           
   J=J+Jwell;
  
   % Add to Jacobian elemenst for derivative of capillary pressure   
   if isfield(fluid, 'pc')
     [dpc, dpc]       = fluid.pc(resSol);
      
      d_p  =      dt.* f_face(:,1).*m_face(:,2).*pcJac;
      
      % d pcflux
      %  --------
      %     dS
      % for a face e_ij gives four contributions to the Jacobi matrix:
      
      dFdSipc1_1 =  d_p.*dpc(ic1)./pv(ic1);  % J_ii
      dFdSipc1_2 =  d_p.*dpc(ic2)./pv(ic1);  % J_ij
      dFdSipc2_1 =  d_p.*dpc(ic1)./pv(ic2);  % J_ji
      dFdSipc2_2 =  d_p.*dpc(ic2)./pv(ic2);  % J_jj      
    
      
   
      J = J + sparse( [ic1;        ic1;        ic2;        ic2], ... %row
                      [ic1;        ic2;        ic1;        ic2], ... %column
            [dFdSipc1_1; -dFdSipc1_2;  -dFdSipc2_1; dFdSipc2_2], ... %value
         G.cells.num, G.cells.num);
      
   end
end
function F = Residual (resSol, resSol_0, dt, fluid, dflux, gflux, pcflux,...
                       findPhaseUpwindCells, G, pv, q,W)
%  DISCRETISATION
%
%  Compute  F given by a residual function Fk for each grid cell
%
%   Fk(s, s0, dt) = s - s0 + dt/pv·[H(s) - max(q,0) - min(q,0)·f],
%
%   H(s) = sum_i [f_i(dflux_i + mo_i*gflux_i)]  (faces i of cell)
%
%  where f_i = mw_i/(mw_i+mo_i), mo_i and mw_i are phase upwind
%  mobilities at face i and dflux and gflux are Darcy flux and gravity flux
%  g(rho1-rho2)*face_normal, respectively. The other flux function f is
%  the fractional flow function evaluated with cell mobilities.
%
%  Advancing an implicit upwind mobility weighted scheme one time step
%  amounts to solving F(s, s0, dt) = 0, given s0 and dt.
%
%  The upwind mobility weighted flux f_i is simply computed by evaluating
%  the water mobility mw in the upwind cell wrt water flux
%
%    wflux = fw(dflux + mo*gflux).
%
%  Similarily, we evaluate mo in the upwind cell wrt oil flux
%
%    oflux = fo(dflux - mw*gflux).
%
%  The upwind cell for each phase at each internal grid face are found
%  using findPhaseUpwindCells, and computeConstData.  The gravity flux is
%  computed in getFlux.

   % Compute cell mobilities
   if isfield(resSol, 'extSat'),
         % Save extremal values of saturation for hysteresis purposes
         % First column contains minima, second maxima
           resSol.extSat(:,1) = min(resSol.s(:,1), resSol_0.extSat(:,1));
           resSol.extSat(:,2) = max(resSol.s(:,1), resSol_0.extSat(:,2));
   end
   mu  = fluid.properties(resSol);
   sat = fluid.saturation(resSol);
   kr  = fluid.relperm(sat, resSol);
   m   = bsxfun(@rdivide, kr, mu);
   
   clear mu sat kr


   internal_faces   = all(G.faces.neighbors~=0, 2);
   ic1  = G.faces.neighbors(internal_faces,1);
   ic2  = G.faces.neighbors(internal_faces,2);
   clear internal_faces

   % Compute water source term
   f_cell = m(:,1) ./ sum(m,2);
   Q      = max(q,0) + min(q,0).*f_cell;
   
   % Calculate source term related to wells. Need to take special care in
   % cases where well has both in and outflux   
   for kk=1:numel(W)
       inperf =  resSol.wellSol(kk).flux<0;
       outperf=  resSol.wellSol(kk).flux>0;
       S_avg  = -resSol.wellSol(kk).flux(inperf).*f_cell(W(kk).cells(inperf));%+W(kk).comp(1);
       influx = -sum(resSol.wellSol(kk).flux(inperf));
       wellhead = sum(resSol.wellSol(kk).flux);
       % compute average saturation in well
       if(wellhead>0)
           S_avg = (sum(S_avg)+wellhead*W(kk).compi(1))./(influx+wellhead);
       else
           S_avg = (sum(S_avg))./(influx);
       end
       Q(W(kk).cells(outperf)) = S_avg*resSol.wellSol(kk).flux(outperf);
       Q(W(kk).cells(inperf))  = resSol.wellSol(kk).flux(inperf).*f_cell(W(kk).cells(inperf));
   end
   clear f_cell

   % Compute the upwind mobility weighted flow of water:
   % For each face k, iw(k) and io(k) is the cell index in which we
   % must evaluate the water and oil mobilities, respectively.
   % always use 3D cells as upwind for 2D/3D bnd
   [iw, io] = findPhaseUpwindCells(m, resSol);
    
   m_face  = [m(iw, 1) m(io,2)];
   m_face = fixCoupledBoundaryMob(G, m_face, resSol, fluid, iw, io);
    
   
   f_face  = m_face(:,1)./sum(m_face,2);
   v_water = f_face.*(dflux+(pcflux(resSol)+ gflux).*(m_face(:,2)));
   clear io iw m m_face f_face

   F   = resSol.s(:,1) - resSol_0.s(:,1);
   F   = F+dt.*(1./pv).*(accumarray([ic1; ic2], [v_water; -v_water]) - Q);

end

function [f, g, pc_flux, pcJac] = getFluxCap(G,state, rock, rho, fluid, ...
      vert_avrg, vert_method, region3D)
% -------------------------------------------------------------------------
%
%                        TECHNICAL DESCRIPTION
%
% -------------------------------------------------------------------------
%
%  g = harmonic average of (n·K·gravity·(rho1-rho2)) in each cell-
%
%  f =
   darcyFaceFlux=state.flux;
   % return total velocity and gravity face flux.
   is_int = all(G.faces.neighbors~=0, 2);
    
   % Internal faces.neighbors
   fn_int = G.faces.neighbors(is_int,:);
   
   % Cell faces
   cf = G.cells.faces(:,1);
   
   [nc, nf] = deal(G.cells.num, sum(double(is_int)));
   cellNo   = rldecode(1 : nc, diff(G.cells.facePos), 2) .';

   % Indices to (internal) half-faces.
   cIntFInx = is_int(G.cells.faces(:,1));

   % Indices of internal face corresponding to cIntFInx.
   globfix  = G.cells.faces(cIntFInx, 1); 

   % rho(1) - rho(2) is always correct (even when rho(2) > rho(1)), because
   % during transport (using, e.g., the 'twophaseUpwBEGrav' transport
   % solver) we only modify the *first* saturation component.
   g_fac   = gravity() * (rho(1) - rho(2));
   dim = size(G.nodes.coords, 2);

   harm_c          = zeros([G.faces.num, 1]);
   harm_g          = zeros([G.faces.num, 1]);
   renum         = zeros([G.faces.num, 1]);
   renum(is_int) = 1 : nf;

   % nKg == n' * K * g *(rho1 - rho2) for all cellfaces.
   [K, r, c] = permTensor(rock, dim);

   assert (size(K,1) == G.cells.num, ...
          ['Permeability must be defined in active cells only.\n', ...
           'Got %d tensors, expected %d (== number of cells).'],   ...
           size(K,1), G.cells.num);
   
   C = zeros(numel(G.faces.neighbors(:,1)), dim);
   
   
   C(is_int,:)= G.cells.centroids(fn_int(:,1),:)...
               -G.cells.centroids(fn_int(:,2),:);
          
   nC = sqrt(sum(C(is_int,:).*C(is_int,:),2));
   
   d_pc = @(cell_pc) (cell_pc(fn_int(:,1),:)- ...
                      cell_pc(fn_int(:,2),:))./nC;

   
   
   nKC = sum(G.faces.normals(cf, r).*K(cellNo,:) .*C(cf,c), 2);
   
   harm_c(is_int) = (2 ./ accumarray(renum(globfix), 1 ./ nKC(cIntFInx)))./nC;

   %bsxfun(@times, K(cellNo,:), C(G.cells.faces(:,1),c), 2)); %hvorfor C?
   % C er på facer.. 
   %harm_c(is_int) = (2 ./ accumarray(renum(globfix), 1 ./ nKC(cIntFInx)))./nC;
   if(~vert_avrg)
      nKg = sum(G.faces.normals(cf, r) .* ...
         bsxfun(@times, K(cellNo,:), g_fac(c)), 2);
      harm_g(is_int) = 2 ./ accumarray(renum(globfix), 1 ./ nKg(cIntFInx));
   else
     % vert_method_3d = 'cells';% not working
     % vert_method_3d = 'cells_pp';% working
     % vert_method_3d = 'topface';
     vert_method_3d = vert_method;
     % remove pc-gravity effects for upper and lower boundary in 2D part of the grid
     if any(region3D)
        cells2D = ~region3D;
        cellF2D = cells2D(cellNo);
     else
        cellF2D = true(size(G.cells.faces,1),1);
     end
     %  hftb = any(bsxfun(@eq,G.cells.faces(:,2),[5,6]),2);
     hftb = any(bsxfun(@eq,G.cells.faces(cellF2D,2),[5,6]),2);
     ftb = G.cells.faces(hftb,1);
     ftb = unique(ftb);
     harm_c( ftb ) = 0;
     switch  vert_method_3d
        case {'cells'}
           % use cellcentroid for defining gravity effert
           % warning('Only simple verion of vertical average is implemented for 3d grid');
           harm_g(is_int) = -norm(g_fac).*harm_c(is_int)...
              .*(G.cells.centroids(G.faces.neighbors(is_int,1),3)...
              -G.cells.centroids(G.faces.neighbors(is_int,2),3))./nC;
           %now gravity flow tru z surfaces
           %harm_g(~is_int)=0;
           harm_g( ftb) = 0;

        case {'cells_pp'}
           % use cellcentroid for defining gravity effert
           nnC=zeros(G.faces.num,1);
           nnC(is_int)=nC;
           hfxy = find(G.cells.faces(:,2) < 5);
           fxy = G.cells.faces(hfxy,1);
           fxy = unique(fxy);
           ifxy = fxy(find(~any(G.faces.neighbors(fxy,:) == 0, 2)));
           xy_faces = ifxy;
           harm_g(xy_faces&faces2D) =  -norm(g_fac).*harm_c(xy_faces).*...
              (G.cells.centroids(G.faces.neighbors(xy_faces,1),3)...
              -G.cells.centroids(G.faces.neighbors(xy_faces,2),3))./nnC(xy_faces);
        case {'topface'}
           % use top surface centroid for defining gravity effect
           % nC and nKC could be calculate from top surface
           % centroid difference in this case
           nnC=zeros(G.faces.num,1);
           nnC(is_int)=nC;
           hfxy = find(G.cells.faces(:,2) < 5);
           fxy = G.cells.faces(hfxy,1);
           fxy = unique(fxy);
           ifxy = fxy(find(~any(G.faces.neighbors(fxy,:) == 0, 2)));
           xy_faces = ifxy;
           cell1=  G.faces.neighbors(xy_faces,1);
           cell2=  G.faces.neighbors(xy_faces,2);
           cellno = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
           cind = find(G.cells.faces(:,2)==5);
           topface(cellno(cind))= G.cells.faces(cind,1);
           harm_g(xy_faces) =  -norm(g_fac).*harm_c(xy_faces).*...
              (G.faces.centroids(topface(cell1),3)...
              -G.faces.centroids(topface(cell2),3))./nnC(xy_faces);
        otherwise
           error('no such vertical average option')
     end
     if any(region3D)
        harm_g3D          = zeros([G.faces.num, 1]);
        nKg = sum(G.faces.normals(G.cells.faces(:,1), r) .* ...
           bsxfun(@times, K(cellNo,:), g_fac(c)), 2);
        harm_g3D(is_int) = 2 ./ accumarray(renum(globfix), 1 ./ nKg(cIntFInx));
        % find 3D faces
        cells3D = region3D;
        cellF3D = cells3D(cellNo);
        faces3D = unique(G.cells.faces(cellF3D, 1));

        %give correct gravity contribution for 3D faces
        harm_g(faces3D) = harm_g3D(faces3D);

        %% TMP: set contribution on 2D/3D boundary to 0
        bndFaces2D3D = G.facesBnd.index;
        harm_c(bndFaces2D3D) = 0;
        harm_g(bndFaces2D3D) = 0;
     end
   end
   
   %nKg=abs(nKg)./sqrt(sum(g_fac.*g_fac,2));
   % Compute harmonic average of nKg(rho1-rho2) on all *internal* faces.
   %harm_c(is_int) = (2 ./ accumarray(renum(globfix), 1 ./ nKC(cIntFInx))).*d_pc./nC;
   %g=(harm_g.*bszfun(@times,g_fac,c)+harm_c*d_pc./nC)./nC;
   pc_flux = @(rSol)  harm_c(is_int).*d_pc(fluid.pc(rSol));
   % precompute some values used to compute dpc/ds in Jacobian
   %pcJac = harm_c(is_int).*G.faces.areas(is_int)./nC;
   pcJac = harm_c(is_int)./nC;
   g = harm_g(is_int);
   f = darcyFaceFlux(is_int);
end

function [f, g] = getFlux(G, rock, rho, darcyFaceFlux)
% -------------------------------------------------------------------------
%
%                        TECHNICAL DESCRIPTION
%
% -------------------------------------------------------------------------
%
%  g = harmonic average of (n·K·gravity·(rho1-rho2)) in each cell-
%
%  f =

   % return total velocity and gravity face flux.
   is_int = all(G.faces.neighbors~=0, 2);

   [nc, nf] = deal(G.cells.num, sum(double(is_int)));
   cellNo   = rldecode(1 : nc, diff(G.cells.facePos), 2) .';

   % Indices to (internal) half-faces.
   cIntFInx = is_int(G.cells.faces(:,1));

   % Indices of internal face corresponding to cIntFInx.
   globfix  = G.cells.faces(cIntFInx, 1);

   % rho(1) - rho(2) is always correct (even when rho(2) > rho(1)), because
   % during transport (using, e.g., the 'twophaseUpwBEGrav' transport
   % solver) we only modify the *first* saturation component.
   g   = gravity() * (rho(1) - rho(2));
   dim = size(G.nodes.coords, 2);

   harm          = zeros([G.faces.num, 1]);
   renum         = zeros([G.faces.num, 1]);
   renum(is_int) = 1 : nf;

   % nKg == n' * K * g *(rho1 - rho2) for all cellfaces.
   [K, r, c] = permTensor(rock, dim);

   assert (size(K,1) == G.cells.num, ...
          ['Permeability must be defined in active cells only.\n', ...
           'Got %d tensors, expected %d (== number of cells).'],   ...
           size(K,1), G.cells.num);

   nKg = sum(G.faces.normals(G.cells.faces(:,1), r) .* ...
             bsxfun(@times, K(cellNo,:), g(c)), 2);

   % Compute harmonic average of nKg(rho1-rho2) on all *internal* faces.
   harm(is_int) = 2 ./ accumarray(renum(globfix), 1 ./ nKg(cIntFInx));

   f = darcyFaceFlux(is_int);
   g = harm(is_int);
end
%--------------------------------------------------------------------------

function constData = computeConstData(G, dflux, gflux, opt)
%
% Given fixed dflux and gflux, we can compute upwind directions for
% phase w given that dflux and gflux have the same sign.  If dflux and
% gflux have opposite signs, then the upwind directions for phase o can be
% computed.
%
% If either dflux or gflux is zero, the phase velocity does not depend on
% the phase mobilities. Then, the sequence induced by the present function
% is as good as any other.
%
% We store these precomputed values to speed up calculations of upwind
% direcions in findFaceMobIx.

   intern     = all(G.faces.neighbors ~= 0, 2);
   v_darcy    = dflux;
   g_vec      = gflux;
   ineighbors = G.faces.neighbors(intern,:);
   constData.ineighbors = ineighbors;

   iw = nan(sum(intern), 1);
   io = iw;

   %Upwind direction for 'water'
   % When v and g have the same sign, the sign of water phase flux is
   % independent of mobilities.
   a = ~(v_darcy < 0) & ~(g_vec < 0);
   b = ~(v_darcy > 0) & ~(g_vec > 0);
   vec = a | b;
   ineighbors       = constData.ineighbors;
   ineighbors(b, :) = ineighbors(b, [2,1]);
   iw(vec) = ineighbors(vec, 1);

   % Upwind dirextion for 'oil'
   % When v and g have the opposite sign, the sign of water phase flux is
   % independent of mobilities.
   a = ~(v_darcy < 0) & ~(g_vec > 0);
   b = ~(v_darcy > 0) & ~(g_vec < 0);
   vec =  a | b;
   ineighbors       = constData.ineighbors;
   ineighbors(b, :) = ineighbors(b, [2,1]);
   io(vec)          = ineighbors(vec, 1);

   constData.v_darcy = v_darcy;
   constData.g_vec   = g_vec;
   constData.io      = io;
   constData.iw      = iw;   
   
   bndCells = reshape(G.faces.neighbors(G.facesBnd.index,:)',[],1);
   
   internNum = int32(intern);   
   internNum(intern>0) = 1:nnz(intern);
   
   constData.bnd3D =  [internNum(G.facesBnd.index) bndCells(opt.region3D(bndCells))];
end

%--------------------------------------------------------------------------

function [IW, IO] = findPhaseUpwindCells(constData, mob)
%Compute upwind cell index for each phase flux (internal to grid).
%
% SYNOPSIS:
%   [iw, io] = findFaceMobIx(constData, mob)
%
%
% PARAMETERS:
%   constData - struct computed in computeConstData
%
%   mob       - Cell mobility. (mob = fluid.mob(resSol)).
%
% RETURNS:
%   iw       - Index for converting from cell mobility to face
%               mobility for water-phase (or primary phase).
%               facemob_w = mob(iw,1).
%
%   io       - Index for converting from cell mobility to face
%               mobility for oil-phase (or secondary phase).
%
% -------------------------------------------------------------------------
%
%                        TECHNICAL DESCRIPTION
%
% -------------------------------------------------------------------------
%
%  For each pair of Darcy flux dflux and gravity flux gflux=g·n(rho1-rho2),
%  we can compute the phase velocities for water (w) and oil (o),
%
%           mw
%     vw = -----  (dflux + mo*gflux)
%          mw+mo
%
%
%           mo
%     vo = -----  (dflux - mw*gflux)
%          mw+mo
%
%
% To use the upwind mobility flux approximation, we evaluate the water
% mobility mw in the upwind direction wrt vw, and the oil mobility mo in
% the upwind direction wrt vo.
%
% This is solved by noting that if dflux and gflux are positive(negative),
% the vw is positive(negative), regardless of the value of mo.  Thus, the
% upwind direction of water is known, mw may be evaluated, vo computed and
% the upwind direction for oil can be found.
%
% If dflux is positive(negative) and gflux is negative(positive), vo is
% positive(negative) and the rest may be computed.

   iw = constData.iw;
   io = constData.io;
   if any(any(isnan([iw, io]))),
      mw = nan(size(iw));
      mo = mw;

      mw(~isnan(iw)) = mob(iw(~isnan(iw)), 1);
      mo(~isnan(io)) = mob(io(~isnan(io)), 2);

      vw = constData.v_darcy +  constData.g_vec.*mo;
      vo = constData.v_darcy -  constData.g_vec.*mw;

      ineighbors = constData.ineighbors;
      ineighbors(vw<0, :) = ineighbors(vw<0, [2,1]);
      iw(isnan(iw)) = ineighbors(isnan(iw), 1);

      ineighbors = constData.ineighbors;
      ineighbors(vo<0, :) = ineighbors(vo<0, [2,1]);
      io(isnan(io)) = ineighbors(isnan(io), 1);
   end   
   
   % always use 3D cell for upwind
%    iw(constData.bnd3D(:,1)) = constData.bnd3D(:,2);
%    io(constData.bnd3D(:,1)) = constData.bnd3D(:,2);
   
   IW = iw;
   IO = io;
end
function q = computeTransportSourceTerm(state, G, wells, src, bc)
%Compute source terms for transport
%
% SYNOPSIS:
%   q = computeTransportSourceTerm(state, G, wells, src, bc)
%
% PARAMETERS:
%   state - Reservoir and well solution structure either properly
%           initialized from functions 'initResSol' and 'initWellSol'
%           respectively, or the results from a call to function
%           'solveIncompFlow'.
%
%   wells - Well structure as defined by function 'addWell'.  May be empty
%           (i.e., W = []) which is interpreted as a model without any
%           wells.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'.  May be empty (i.e., src = []) which is
%           interpreted as a reservoir model without explicit sources.
%
%   bc    - Boundary condition structure as defined by function 'addBC'.
%           This structure accounts for all external boundary conditions
%           to the reservoir flow.  May be empty (i.e., bc = []) which is
%           interpreted as all external no-flow (homogeneous Neumann)
%           conditions.
%
% RETURNS:
%   q     - Aggregate source term contributions (per grid cell) suitable
%           for passing to a transport solver.  Measured in units of m^3/s.
%
% SEE ALSO:
%   `twophaseJacobian`, `computeTimeOfFlight`.

% TODO:
%   - implement gravity effects for pressure boundary and wells

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

% $Date: 2012-10-05 10:15:39 +0200 (Fri, 05 Oct 2012) $
% $Revision: 10006 $

   qi = [];  % Cells to which sources are connected
   qs = [];  % Actual strength of source term (in m^3/s).

   if ~isempty(wells),
      [i, s] = contrib_wells(wells, state.wellSol);
      qi = [qi; i];
      qs = [qs; s];
   end

   if ~isempty(src), assert (~isempty(src.sat))
      [i, s] = contrib_src(src);
      qi = [qi; i];
      qs = [qs; s];
   end
   
   if ~isempty(bc), assert (~isempty(bc.sat))
      is_int = all(double(G.faces.neighbors) > 0, 2);
      [i, s] = contrib_bc(G, state, bc, is_int);
      qi = [qi; i];
      qs = [qs; s];
   end

   %-----------------------------------------------------------------------
   % Assemble all source and sink contributions to each affected cell. ----
   %
   q = sparse(qi, 1, qs, G.cells.num, 1);
end

%--------------------------------------------------------------------------

function [qi, qs] = contrib_wells(W, wellSol)
   % Wells as defined by 'addWell'.

   nperf = cellfun(@numel, { W.cells }) .';

   qi = vertcat(W.cells);
   qs = vertcat(wellSol.flux);

   % Injection perforations have positive flux (into reservoir).
   %
   comp      = rldecode(vertcat(W.compi), nperf);
   inj_p     = qs > 0;
   qs(inj_p) = qs(inj_p) .* comp(inj_p,1);
end

%--------------------------------------------------------------------------

function [qi, qs] = contrib_src(src)
   % Explicit sources defined by (e.g.) 'addSource'.

   qi = src.cell;
   qs = src.rate;

   % Injection sources have positive rate into reservoir.
   %
   in = find(src.rate > 0);
   if ~isempty(in),
      qs(in) = qs(in) .* src.sat(in,1);
   end
end

%--------------------------------------------------------------------------

function [qi, qs] = contrib_bc(G, state, bc, is_int)
   % Contributions from boundary conditions as defined by 'addBC'.

   qs = zeros([G.faces.num, 1]);
   dF = false([G.faces.num, 1]);

   isDir = strcmp('pressure', bc.type);
   isNeu = strcmp('flux',     bc.type);

   dF(bc.face(isDir))      = true;
   cfIx                    = dF(G.cells.faces(:,1));

   cflux = faceFlux2cellFlux(G, state.flux);
   qs(G.cells.faces(cfIx,1)) = -cflux(cfIx);
   qs(bc.face(isNeu))      = bc.value(isNeu);

   % Injection BC's have positive rate (flux) into reservoir.
   %
   is_inj = qs > 0;
   if any(is_inj),
      qs(is_inj) = qs(is_inj) .* bc.sat(is_inj(bc.face), 1);
   end

   is_outer = ~is_int;

   qi = sum(G.faces.neighbors(is_outer,:), 2);
   qs = qs(is_outer);
end

function varargout = fixCoupledBoundaryMob(G, m_face, state, fluid, iw, io, varargin)
    % The upstream mobility must be altered over the edge of the special
    % domains.
    fint = all(G.faces.neighbors,2);
    loc = nan(G.faces.num,1);
    loc(fint) = 1:sum(fint);
    % Look up local indices
    loc_bnd2d = loc(G.facesBnd.index);
    % Find upstream cells on the 2D part of the boundary fort both phases
    w2d = ismember(iw(loc_bnd2d), G.facesBnd.cells2D);
    o2d = ismember(io(loc_bnd2d), G.facesBnd.cells2D);
    
    if ~any(w2d | o2d)
        return;
    end
    
    % Find the cells which are neighbors to a 3D region (i.e. both 2D
    % and 3D faces)

    bndCells2D = G.facesBnd.cells2D;

    % Number of fine cells in the full grid in these columns
    numSubCells = (G.cells.columnPos(G.cells.mapTopSurface(bndCells2D)+1)-...
        G.cells.columnPos(G.cells.mapTopSurface(bndCells2D)));

    % Find 2D height values via the topsurfacegrid
    [s_full h_top] = normalizeValuesVE(G.topSurfaceGrid, state, fluid, 'CoupledGrid', G);
    % Once we have height, we can find filldegree of the cells on the
    % boundary.
    [n, t] = fillDegree(h_top, G.topSurfaceGrid);
    % Fill degree of the partially filled cells on the edge of the 3/2D
    % split.
    fill2d = t(G.facesBnd.cells2D);

    n2D = numel(G.facesBnd.cells2D);
    Gt = G.topSurfaceGrid;
    cP = Gt.cells.columnPos;
    partial = false(numel(w2d),1);
    filled  = partial;
    % Find partial, filled and empty cells.
    for i = 1:n2D
        ind = G.facesBnd.cells2D(i);
        cells = Gt.columns.cells(cP(ind):cP(ind+1)-1);
        partcell = G.facesBnd.map3D == cells(n(ind)+1);
        partial(partcell) = true;
        filled(sum(numSubCells(1:i-1)) + 1: find(partcell)) = true;
    end
    % Empty is true where there is no CO2
    % Filled is true where s_CO2 == 1
    % Partial is true where s_CO2 is not 0 or 1
    empty = ~filled;
    filled = filled & ~partial;
    fill = rldecode(fill2d, numSubCells);


    mu = fluid.properties(state);
    
    % This area is only filled with CO2 and will therefore have relperm = 1
    % for CO2 and 0 for other phase.
    m_face(loc_bnd2d(w2d & filled),1) = 1/mu(1);
    m_face(loc_bnd2d(o2d & filled),2) = 0;
    
    % The second phase is the opposite case
    m_face(loc_bnd2d(w2d & empty),1) = 0;
    m_face(loc_bnd2d(o2d & empty),2) = 1/mu(2);
    
    % For the partial saturation, we use a linear combination of the fill
    % degree (which is already \in [0,1]...)
    m_face(loc_bnd2d(w2d & partial),1) =    fill(w2d & partial )./mu(1);
    m_face(loc_bnd2d(o2d & partial),2) = (1-fill(o2d & partial))./mu(2);

    varargout{1} = m_face;
    if nargout > 1
        % d(mob)/dx is only nonzero in the partial cell since this is the
        % only region sensitive to changes in saturation (the cells above
        % and below are empty and full as long as the cells are not
        % intersected by the plume boundary).
        %
        % The partially filled cells will have d(mob_i)/ds =
        % d(mob_i)/d(filldegree) = 1/mu_co_2 and -1/mu_w
        dm_face = varargin{1};
        dm_face(loc_bnd2d(w2d & filled),1) = 0;
        dm_face(loc_bnd2d(w2d & empty),1) = 0;
        dm_face(loc_bnd2d(w2d & partial),1) = 1/mu(1);

        dm_face(loc_bnd2d(o2d & filled),2) = 0;
        dm_face(loc_bnd2d(o2d & empty),2) = 0;
        dm_face(loc_bnd2d(o2d & partial),2) = -1/mu(2);

        varargout{2} = dm_face; 
    end
end
