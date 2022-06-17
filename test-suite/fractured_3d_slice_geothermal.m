function setup = fractured_3d_slice_geothermal(varargin)
% 

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
    % One-line description
    description = '';
    
    % Optional input arguments
    options = struct('cartDims', [50,1,50], ... % Number of cells in x- and y-directions
                     'physDims', [2,1,20] , ...
                     'numFrac' , 20       , ...
                     'apertureFrac', 1*milli*meter, ...
                     'time'    , 14*day   , ... % Total injection time
                     'dt'      , 6*hour   );    % Brooks-Corey relperm exponent
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props geothermal compositional

    gravity reset on
    
    % Cartesian grid
    G    = computeGeometry(cartGrid(options.cartDims, options.physDims));
    rock = makeRock(G, 1*milli*darcy, 0.05);
    rock = addThermalRockProps(rock, 'lambdaR', 2*watt/(meter*Kelvin)       , ...
                                     'rhoR'   , 2700*kilogram/meter^3       , ...
                                     'CpR'    , 1000*joule/(kilogram*Kelvin));
    K0   = 273.15*Kelvin;
    % Make fluid
    fluid = initSimpleADIFluid(                     ...
                   'phases', 'W'                  , ... % Water only
                   'mu'    , 0.5*centi*poise      , ... % Viscosity
                   'rho'   , 1000*kilogram/meter^3, ... % Reference density
                   'c'     , 4.4e-10/Pascal       , ... % Compressibility
                   'pRef'  , 1*atm                );    % Ref pressure
    % Add thermal properties
    fluid = addThermalFluidProps(fluid            , ...
                'Cp'     , 4.2*joule/(gram*Kelvin), ... % Heat capacity
                'lambdaF', 0.6*watt/(meter*Kelvin), ... % Thermal cond
                'cT'     , 207e-6/Kelvin          , ... % Thermal expansion
                'TRef'   , K0 + 10*Kelvin         );    % Reference temp
    % Make model
    [ii, jj, kk] = gridLogicalIndices(G);
    
    cells = any(kk == (5:5:numel(kk)-5), 2);
    rock.perm(cells) = options.apertureFrac^2/12;
    rock.poro(cells) = 0.5;
    
    rock.perm = rock.perm*0 + options.apertureFrac^2/12*0.1;
    rock.poro = rock.poro*0 + 0.5;
    
    % Construct two-phase model with oil and water
    model = GeothermalModel(G, rock, fluid);
    
    model = model.validateModel();
    model.FacilityModel.implicitConnectionDP = true;
    
    % Wells
    rate = 0.5*litre/second;
    bhp  = 50*barsa;
    W = [];
    W = verticalWell(W, G, rock, 1, 1, [], 'type', 'rate', 'val', rate, 'comp_i', [1,0]);
    W = verticalWell(W, G, rock, options.cartDims(1), 1, [], 'type', 'bhp', 'val', bhp, 'comp_i', [1,0]);
    W = addThermalWellProps(W, G, rock, fluid, 'T', convertFromCelcius(95));
    
%     for i = 1:numel(W)
%         W(i) = convert2MSWell(W(i), 'G', G);
%     end
    
    % Schedule
    schedule = simpleSchedule(rampupTimesteps(options.time, options.dt), 'W', W);
    
    % Initial state
    state0 = initResSol(G, 1,1);
    [state0.pressure, state0.T] = initializeGeothermalEquilibrium(model);
    % Pack setup
    
    % Plotting
    plotOptions = {'View', [0,0], ...
                   'Projection', 'orthographic', ...
                   'Size', [700, 500], ...
                   'PlotBoxAspectRatio', [1,1,3]};
    
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end