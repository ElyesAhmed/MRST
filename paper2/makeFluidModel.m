function f  = makeFluidModel(aquifer, varargin)
% Make an ADI fluid for the 1D antiform aquifer model
%
% SYNOPSIS
%   f = makeFluidModel(aquifer)
%   f = makeFluidModel(aquifer, 'pn1','pv1')
%
% PARAMETERS:
%    - aquifer     : structure representing the aquifer.  Has two fields:
%                    - 'G'      - the original 3D grid
%                    - 'Gt'     - a top surface grid representation of 'G'
%                    - 'rock'   - a structure containing the rock parameters
%                                 (permeability and porosity) for the original
%                                 3D grid
%                    - 'rock2D' - a structure containing the rock parameters
%                                 for the vertically averaged grid.
%                    - 'W'      - wells structure
%
%   'pn'/pv - List of property names/property values that decide which type
%   of model to use. These are the parameters and values:
%    - residual    : true/false. Default: true
%    - dissolution : true/false. Default: false
%    - fluidType   : 'simple', 'integrated', 'sharp interface',
%                    'linear cap', 'S table', 'P-scaled table',
%                    'P-K-scaled table'. Default: 'sharp interface'
%   For sharp-interface models, there are two additional parameters:
%    - top_trap    : vector describing average height of sub-scale
%                    topographic traps. Default: []
%    - surf_topo   : string describing the type of sub-scale topography:
%                    no variation ('smooth'), accretion layer model ('inf
%                    rough'), sinusoidal geometry ('sinus'), step function
%                    ('square'). Default: 'smooth'
%    - only_pvt    : do not give relperm functions, useful for plotting
%                    values are true or false
%    - co2_type    : choose between different pvt calculations, all options
%                    else than default need third party pvt calculations
%                    and mrst-external loaded
%
% SEE ALSO:
%   makeAquiferModel

opt = struct(...
    'residual'    , true              , ...
    'dissolution' , false             , ...
    'fluidType'   , 'sharp interface' , ...
    'top_trap'    , []                , ...
    'surf_topo'   , 'smooth'          , ...
    'temp_grad'   , 30 / (kilo*meter) , ...
    'surf_temp'   , 12                , ... % in Celsius
    'only_pvt'    , false             , ...% used for plotting
    'co2_type'    , 'with_phase_boundary' ...%
    );
opt = merge_options(opt,varargin{:});

%% Set parameters
mu = [6e-5, 8e-4] * Pascal * second; 
rho = [760 1100] * kilogram / meter^3; 
if (opt.residual)
   [res_gas, res_water] = deal(.21,.11); 
else
   [res_gas, res_water] = deal(0); 
end

f = initSimpleADIFluid('mu'  , mu([2 2 1])  ,...
                       'rho' , rho([2 2 1]) ,...
                       'n'   , [1 1 1]); 

for wfield = {'krO', 'krW', 'krG', 'pcWG', 'pcOG'}
   if (isfield(f, wfield))
      f = rmfield(f, wfield);
    end
end

f.pvMultR = @(p) 1 + (1e-5 / barsa) * (p - 100 * barsa); 
f.bW      = @(p, varargin) 1 + (4.3e-5 / barsa) * (p - 100 * barsa); 
f.BW      = @(p, varargin) 1 ./ f.bW(p); 
T_res     = aquifer.Gt.cells.z * opt.temp_grad + (274 + opt.surf_temp); % tempratures

switch opt.co2_type
  case 'with_phase_boundary'
    f.bG       = boCO2(T_res, f.rhoGS); 
  case 'coolprops'
    co2        = coolProps2Pvt('co2'); 
    f.bG       = @(p) co2.density(p, T_res) ./ f.rhoGS; 
  case 'coolprops_table'
    co2        = coolProps2Pvt('co2'); 
    co2        = pvt2Table({co2}, linspace(100, 500, 40) * barsa, ...
                    linspace(min(T_res) - 10, max(T_res) + 10, 40)); 
    co2        = co2{1}; 
    f.bG       = @(p) co2.density(p, T_res) ./ f.rhoGS; 
  case  'coolprops_linear' 
    assert(norm(gravity)>0)
    p_ref      = aquifer.Gt.cells.z * f.rhoWS * norm(gravity); 
    co2        = coolProps2Pvt('co2'); 
    [p_ref, i] = max(p_ref); 
    co2        = pvt2Linear(co2, p_ref, T_res(i)); 
    f.bG       = @(p) co2.density(p, T_res) ./ f.rhoGS; 
    otherwise
        error('no such co2 pvt')
end

f.BG      = @(p) 1./f.bG(p);
f.surface_tension = 30e-3;

% Constants and definitions shared by several models
drho    = f.rhoWS - f.rhoGS; 
g       = norm(gravity); 
C       = 0.4 * g * max(aquifer.Gt.cells.H) * drho; 
alpha   = 0.5; 
beta    = 3.0; 
samples = 100; 
tabSw   = linspace(0, 1, 10)'; 
tabW    = struct('S', 1 - tabSw, 'kr', tabSw, 'h', []); 

if opt.dissolution
    f.dis_rate = 5e-11;
    f.dis_max  = 0.03;
    f.muW      = @(po,rs,flag,varargin) f.muW(po);
    f.rsSat    = @(po,rs,flag,varargin) (po*0+1)*f.dis_max;
end

if(~opt.only_pvt)
   switch opt.fluidType
     case 'simple'
       f.krG     = @(sg, varargin) sg; 
       f.krW    = @(sw, varargin) sw; 
       f.pcWG    = @(sg, p, varargin) ...
                    g * ( f.rhoWS .* f.bW(p) - f.rhoGS .* f.bG(p) ) .* ...
                    (sg) .* aquifer.Gt.cells.H; 
       f.res_gas = 0; 
       f.res_water = 0; 
       f.invPc3D = @(p) 1 - (sign(p + eps) + 1) / 2; 
       f.kr3D    = @(s) s; 
       
     case 'integrated'
       f = addVERelpermIntegratedFluid(f ,...
           'res_water'   , res_water     ,...
           'res_gas'     , res_gas       ,...
           'Gt'          , aquifer.Gt    ,...
           'kr_pressure' , true          ,...
           'int_poro'    , false         ,...
           'rock'        , rock); 
       
     case 'sharp interface'
       f = addVERelperm(f ,...
           aquifer.Gt                 ,...
           'res_water'   , res_water      ,...
           'res_gas'   , res_gas      ,...
           'top_trap'  , opt.top_trap ,...
           'surf_topo' , opt.surf_topo); 
       
     case 'linear cap'
       f = addVERelpermCapLinear(f, 0.2 * g * max(aquifer.Gt.cells.H) * drho, ...
           'res_gas'     , res_gas            ,...
           'res_water'   , res_water          ,...
           'beta'        , 2                  ,...
           'H'           , aquifer.Gt.cells.H ,...
           'kr_pressure' , true); 

     case 'S table'
       table_co2_1d = makeVEtables(...
           'invPc3D'    , @(p) max((C ./ (p + C)).^(1 / alpha) , res_water) ,...
           'is_kscaled' , false                                           ,...
           'kr3D'       , @(s) s.^beta                                    ,...
           'drho'       , drho                                            ,...
           'Gt'         , aquifer.Gt                                      ,...
           'samples'    , samples); 
       f = addVERelperm1DTables(f                                         ,...
           'res_water'   , res_water                                      ,...
           'res_gas'     , res_gas                                        ,...
           'height'      , aquifer.Gt.cells.H                             ,...
           'table_co2'   , table_co2_1d                                   ,...
           'table_water' , tabW); 
       
     case 'P-scaled table'
       table_co2_1d = makeVEtables(...
           'invPc3D'    , @(p) max((C ./ (p + C)).^(1 / alpha) , res_water) ,...
           'is_kscaled' , false                                           ,...
           'kr3D'       , @(s) s.^beta                                    ,...
           'drho'       , drho                                            ,...
           'Gt'         , aquifer.Gt                                      ,...
           'samples'    , samples); 
       f = addVERelperm1DTablesPressure(f                                 ,...
           'res_water'     , res_water                                        ,...
           'res_gas'     , res_gas                                        ,...
           'height'      , aquifer.Gt.cells.H                             ,...
           'table_co2'   , table_co2_1d                                   ,...
           'table_water' , tabW                                           ,...
           'kr_pressure' , true); 
       
     case 'P-K-scaled table'
       C = 1; 
       kscale = sqrt(aquifer.rock2D.poro ./ (aquifer.rock2D.perm)) * f.surface_tension; 
       table_co2_1d = makeVEtables(...
           'invPc3D'    , @(p) max((C ./ (p + C)).^(1 / alpha) , res_water) ,...
           'is_kscaled' , true                                            ,...
           'kr3D'       , @(s) s.^beta                                    ,...
           'drho'       , drho                                            ,...
           'Gt'         , aquifer.Gt                                      ,...
           'samples'    , samples                                         ,...
           'kscale'     , kscale); 
       f = addVERelperm1DTablesPressure(f     ,...
           'res_water'     , res_water            ,...
           'res_gas'     , res_gas            ,...
           'height'      , aquifer.Gt.cells.H ,...
           'table_co2'   , table_co2_1d       ,...
           'table_water' , tabW               ,...
           'rock'        , aquifer.rock       ,...
           'kr_pressure' , true); 
       
     otherwise
       error('No such fluid case')
   end
else
   f.res_gas = res_gas; 
   f.res_water = res_water; 
end
end
