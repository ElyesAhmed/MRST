function fluidPlotPanelAD(model, varargin)
% Create a interactive plotting panel for a given model that shows
% different fluids properties.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
    opt = struct('pressureRange', []);
    opt = merge_options(opt, varargin{:});

    if isempty(opt.pressureRange)
        p0 = max(model.minimumPressure, 0);
        p1 = min(model.maximumPressure, 1000*barsa);
        opt.pressureRange = subdiv(p0, p1);
    end
    % Make a figure that's wider than the default.
    df = get(0, 'DefaultFigurePosition');
    fh = figure('Position', df.*[1 1 1.75 1]);
    % We want to ensure that the lines are nice and pretty even if the user
    % has messed with the defaults.
    set(fh, 'Renderer', 'painters');
    
    % Options can alter the amount of space the plot itself takes up.
    lm = 0.1;
    pw = .7;
    % Somewhat magic numbers because the gui in matlab has some magic
    % constants itself.
    plotaxis  = subplot('Position', [.75*lm, lm, pw, 1-2*lm]);
    ctrlpanel = uipanel('Position', [lm+pw, .25*lm, 1-1.25*lm-pw,  1-1.25*lm], ...
                        'Title',  'Property');
    true3ph = model.water && model.oil && model.gas;
    names = {'Viscosity',...
             'Relative permeability', ...
             'Rock compressibility', ...
             'Densities', ...
             'Capillary pressure', ...
             'Max Rs/Rv', ...
             '3ph-relperm: Water', ...
             '3ph-relperm: Oil', ...
             '3ph-relperm: Gas' ...
             };
    functions = {@(name) plotViscosity(model, name), ...
                 @(name) plotRelperm(model, name),...
                 @(name) plotPVMult(model, name),...
                 @(name) plotDensity(model, name), ...
                 @(name) plotCapillaryPressure(model, name), ...
                 @(name) plotMaxR(model, name), ...
                 @(name) plot3phRelPerm(model, name, 1), ...
                 @(name) plot3phRelPerm(model, name, 2), ...
                 @(name) plot3phRelPerm(model, name, 3) ...
                };
    
    active = [true; ...
              true; ...
              isfield(model.fluid, 'pvMult'); ...
              true; ...
              isfield(model.fluid, 'pcOW') || isfield(model.fluid, 'pcOG'); ...
              model.disgas || model.vapoil; ...
              true3ph;
              true3ph;
              true3ph];
              
              
    names = names(active);
    functions = functions(active);
    
    phases = {'Water', 'Oil', 'Gas'};
    for it = 1:numel(phases)
        ph = phases{it};
        if model.(lower(ph));
            names = [names, [ph, ' viscosity'], ...
                            [ph, ' b-factor']];
            functions = [functions, {@(name) plotStuff(model, {['mu', ph(1)]}, name), ...
                                     @(name) plotStuff(model, {['b', ph(1)]}, name)}];
        end
    end
    
    propsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'listbox',...
              'String', names, 'Callback', @drawPlot, ...
              'Position',[0 0 1 1]);
    
          

    function drawPlot(src, event)
        axis(plotaxis);
        colorbar off;
        axis normal
        view(0, 90);
        ix = get(propsel, 'Value');
        fn = functions{ix};
        name = names{ix};
        
        fn(name);
    end

    drawPlot([], []);


    function plotStuff(model, fields, plottitle)
        f = model.fluid;
        p = opt.pressureRange;
        s = subdiv(0, 1);

        rsMax = 0;
        rvMax = 0;
        if model.disgas
            rsMax = f.rsSat(p);
        end
        if model.vapoil
            rvMax = f.rvSat(p);
        end
        n = sum(model.getActivePhases);

        legflag = false(size(fields));
        legh = zeros(size(fields));
        ctr = 0;
        yl = '';
        cla;
        hold on
        nf = numel(fields);
        colors = lines(nf);
        for i = 1:nf
            fn = fields{i};
            legflag(i) = true;
            ctr = ctr + 1;
            switch(lower(fn))
                case {'krw', 'krg', 'krog', 'kro', 'krow'}
                    data = f.(fn)(s);
                    x = s;
                    xl = 'Saturation';
                case {'pcow', 'pcog'}
                    data = f.(fn)(s);
                    x = s;
                    xl = 'Saturation';
                case {'bo', 'bg', 'bw'}
                    [x, data, ok] = evalSat(model, f, fn, p, rsMax, rvMax);
                    x = x/barsa;
                    xl = 'Pressure (barsa)';
                case {'rhoo', 'rhog', 'rhow'}
                    bsub = ['b', fn(end)];
                    rho = f.(['rho', fn(end), 'S']);
                    [x, b, ok] = evalSat(model, f, bsub, p, rsMax, rvMax);
                   
                    data = b*rho;
                    x = x/barsa;
                    xl = 'Pressure (barsa)';
                    yl = 'Density [kg/m^3]';
                case {'muw', 'muo', 'mug'}
                    [x, data, ok] = evalSat(model, f, fn, p, rsMax, rvMax);
                    x = x/barsa;
                    xl = 'Pressure (barsa)';
                case {'rssat', 'rvsat'}
                    data = f.(fn)(p);
                    x = p/barsa;
                    xl = 'Pressure (barsa)';
                case {'pvmultr'}
                    data = f.(fn)(p);
                    x = p/barsa;
                    xl = 'Pressure (barsa)';
            end
            if size(data, 2) > 1
                for j = 1:size(data, 2)
                    o = ok(:, j);
                    % Ok are saturated lines, draw as thick lines
                    h = plot(x(o, j), data(o, j), 'linewidth', 2, 'color', colors(i, :));
                    % Not ok are undersaturated, draw as thin lines
                    plot(x(~o, j), data(~o, j), '--', 'linewidth', 1, 'color', colors(i, :));
                end
            else
                h = plot(x, data, 'linewidth', 2, 'color', colors(i, :));
            end
            legh(i) = h(1);
        end
        
        grid on
        legend(legh, fields(legflag))
        xlabel(xl)
        ylabel(yl);
        if nargin == 3
            title(plottitle)
        end
    end

    function plotViscosity(model, name)
        plotStuff(model, {'muW', 'muO', 'muG'});
        title('Viscosity');
    end

    function plotRelperm(model, name)
        krnames = {};
        if model.water
            krnames = [krnames, 'krW'];
        end
        if model.oil
            if model.water && model.gas
                krnames = [krnames, 'krOW', 'krOG'];
            else
                krnames = [krnames, 'krO'];
            end
        end
        if model.gas
            krnames = [krnames, 'krG'];
        end
        plotStuff(model, krnames);
        title(name);
    end

    function plotDensity(model, name)
        plotStuff(model, {'rhoW', 'rhoO', 'rhoG'});
        title(name);
    end

    function plotCapillaryPressure(model, name)
        cnames = {};
        fld = {'pcOW', 'pcOG'};
        for i = 1:numel(fld)
            if isfield(model.fluid, fld{i})
                cnames = [cnames, fld{i}];
            end
        end

        plotStuff(model, cnames);
        title('Capillary pressure');
    end

    function plotMaxR(model, name)
        rnames = {};
        if model.disgas
            rnames = [rnames, 'rsSat'];
        end
        if model.vapoil
            rnames = [rnames, 'rvSat'];
        end
        plotStuff(model, rnames);
        title(name);
    end

    function plotPVMult(model, name)
        plotStuff(model, {'pvMultR'});
        title(name);
    end
    
    function plot3phRelPerm(model, name, ix)
        cla;
        s = subdiv(0, 1, 50);
        [x, y] = meshgrid(s);
        [krW, krO, krG] = deal(zeros(size(x)));
        
        for i = 1:size(x, 1)
            xi = x(i, :);
            yi = y(i, :);
            [krW(i, :), krO(i, :), krG(i, :)] = model.relPermWOG(xi, 1 - xi - yi, yi, model.fluid);
        end
        
        if ix == 1
            surf(x, y, krW)
            title('Water relative permeability')
        elseif ix == 2
            surf(x, y, krO)
            title('Oil relative permeability')
        else
            surf(x, y, krG)
            title('Gas relative permeability')
        end
        caxis([0, 1]);
        shading interp
        axis equal tight
        % view(0, 90)
        view(-35, 45);
        xlabel('S_w')
        ylabel('S_g')
        colorbar
        legend off
    end

    function x = subdiv(start, stop, n)
        if nargin < 3
            n = 100;
        end
        dx = (stop - start)/n;
        x = (start:dx:stop)';
    end

    function [x, y, ok] = evalSat(model, f, fn, x, rsMax, rvMax)
        ok = true(size(x));
        if checkBO(model)
            if any(strcmpi(fn, {'muo', 'bo'})) && model.disgas
                mrs = max(rsMax);
                rs = 0:mrs/10:mrs;
                [x, rs_g] = meshgrid(x, rs);
                rssat = zeros(size(x));
                for i = 1:size(x, 1)
                    rssat(i, :) = f.rsSat(x(i, :));
                end

                saturated = rs_g >= rssat;
                rs_g(saturated) = rssat(saturated);
                ok = saturated';
                y = f.(fn)(x, rs_g, saturated)';
                x = x';
            elseif any(strcmpi(fn, {'mug', 'bg'})) && model.vapoil
                mrs = max(rvMax);
                rv = 0:mrs/10:mrs;
                [x, rs_g] = meshgrid(x, rv);
                rvsat = zeros(size(x));
                for i = 1:size(x, 1)
                    rvsat(i, :) = f.rvSat(x(i, :));
                end

                saturated = rs_g >= rvsat;
                rs_g(saturated) = rvsat(saturated);
                ok = saturated';
                y = f.(fn)(x, rs_g, saturated)';
                x = x';
            else
                y = f.(fn)(x);
            end
        else
            y = f.(fn)(x);
        end
    end

    function ind = checkBO(model)
        ind = isa(model, 'ThreePhaseBlackOilModel') &&...
               (model.disgas || model.vapoil);
    end
end