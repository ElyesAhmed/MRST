classdef postProcessDiagnostics < handle
    properties
        Figure
        Axes3D
        Axes2DL
        Axes2DR
        colorBar
        colorHAx
        G
        Gs
        WellPlot
        Data
        Menu
        Props
        Measures
        Allocation
        Distributions
        Patch
        maxTOF
        currentInteractionRegion = []
        currentFilterRegion      = []
        currentDiagnostics
        interactiveModes
        camlight
        outlineGrid
        info = []
        layoutParams = struct('menuWidth',        300, ...
                              'distAxisHeight',   150, ...
                              'itemSpace',         25, ...
                              'includeToolTips', true);
    end
    
    methods
        function d = postProcessDiagnostics(varargin)
            opt = struct('style', 'default', 'steps', [], 'maxTOF', 500*year);
            if mod(numel(varargin), 2) == 1  % file-name provided
                filenm = varargin{1};
                [pth, nm, ext] = fileparts(filenm);
                if ~strcmp(ext, 'EGRID')
                    filenm = fullfile(pth, [nm, '.', 'EGRID']);
                end
                opt = merge_options(opt, varargin{2:end});
            else
                [fn, pth] = uigetfile('*.EGRID', 'Select ECLIPSE output grid file (EGRID)');
                if ~isempty(fn) && ischar(fn)
                    filenm = fullfile(pth, fn);
                else
                    return;
                end
                opt = merge_options(opt, varargin{:});
            end
            assert(exist(filenm, 'file')>0, sprintf('Unable to find file %s', filenm));
                
            % Select which time-steps to include
            [pth, nm] = fileparts(filenm);
            casenm = fullfile(pth, nm);
            rsspec = processEclipseRestartSpec(casenm, 'all');
            if ~isempty(opt.steps)
                steps = opt.steps;
            else
                steps = uiPreSelectTimeSteps(rsspec);
            end
       
            % Setup data for selected steps and compute diagnostics
            d.maxTOF = opt.maxTOF;
            [d.G, d.Data, d.Gs, valid_ix] = readAndPrepareForPostProcessor(casenm, steps, rsspec);
            d.Data = computeDiagnostics(d.Gs, d.Data, d.maxTOF, []);
            try
                d.Data.summary = readEclipseSummaryUnFmt(casenm);
            catch
                d.Data.summary = [];
                warning('Not able to read summary for selected case ...\n')
            end
            
            % ------ Set unique colors for wells --------------------------
            % Avoid first color (black) and add an extra color representing
            % potential contributions from the reservoir
            [nInj,nProd] = size(d.Data.wellComunication);
            cmap = tatarizeMap(nInj+nProd+2);
            d.Data.injColors  = cmap([2:nInj+1 end],:);
            d.Data.prodColors = cmap(nInj+2:end,:);
            
            % ------ Create figure window ---------------------------------
            screensize = get( groot, 'Screensize' );
            wsize = .75*screensize(3:4);
            wsize = [screensize(3:4)-wsize*1.1, wsize];
            d.Figure = limitedToolbarFigure('Position', wsize);
            
            % ------ Selectors for time steps and wells -------------------
            itemOpts = {'Parent', d.Figure, 'Visible','off', 'style', opt.style};
            
            tStepStr = getTStepStrings(d.Data.time.cur, d.Data.time.prev);
            selector3D.tsel = timeStepSelector('tSteps', tStepStr, itemOpts{:});
            selector3D.wsel = wellSelector(...
                    'injectors', {d.Data.diagnostics(1).WP.inj.name}, ...
                    'producers', {d.Data.diagnostics(1).WP.prod.name}, itemOpts{:});
            selector3D.wsel.communicationMatrix = d.Data.wellComunication;
            
            % ------ Selector for properties in 3D plot -------------------
            d.Props = struct('static',  struct('name', {{d.Data.static.name}},  'limits', {{d.Data.static.limits}}), ...
                             'dynamic', struct('name', {{d.Data.dynamic.name}}, 'limits', {{d.Data.dynamic.limits}}), ...
                             'diagnostics', struct('name', {{'TOF forward', 'TOF backward', ...
                                                             'Residence time', 'Tracer forward', ... 
                                                             'Tracer backward', 'Tracer product', ...
                                                             'Sweep regions', 'Drainage regions', ...
                                                             'First arrival forward', 'First arrival backward'}}, ...
                                                   'limits', {{[0 1.001*d.maxTOF/year], [0 1.001*d.maxTOF/year], ...
                                                               [0 2.001*d.maxTOF/year], [0 1], [0 1], [0 1], [], [], ...
                                                               [0 1.001*d.maxTOF/year], [0 1.001*d.maxTOF/year]}}), ...
                             'computed', struct('name', {{}}, 'limits', {{}}));
            d.currentDiagnostics = emptyDiagnostics(d);
            selector3D.psel = propertyDisplaySelector('props', d.Props, itemOpts{:}, 'includeLogSwitch', true);
            
            % ------ Selector for property filter in 3D -------------------
            selector3D.fsel = propertyDisplaySelector( ...
                                'Title', 'Property filter', 'includeFilter', true, ...
                                'props', d.Props, 'includePlayer', true, ...
                                itemOpts{:},...
                                'includeEnableSwitch', true, 'includeLogSwitch', true, ... 
                                'startPlayCallback', @d.startPlayCallback, ...
                                'stopPlayCallback', @d.stopPlayCallback);
            selector3D.fsel.Min = d.Data.static(1).limits(1);
            selector3D.fsel.Max = d.Data.static(1).limits(2);
            
            % ------ Selector for heterogeneity measures ------------------
            d.Measures = {{'none','Well connections', 'F-Phi plot', ...
                'Sweep efficiency', 'Lorenz coefficient'}};
            selector2D.msel = dynamicMeasureSelector('props', d.Measures, itemOpts{:});
            
            % ------ Selector for well allocation -------------------------
            d.Allocation = {{'none','Injector volumes', ...
                'Injector allocation', 'Injector profile', ...
                'Producer volumes', 'Producer allocation', 'Producer profile'}};
            selector2D.asel = dynamicMeasureSelector(...
                'Title','Well allocations', ...
                'props', d.Allocation, ...
                'includeAvgSwitch', true, itemOpts{:});
            
            % ------ Selector for summary output  -------------------------
            selector2D.ssel = summarySelector(d.Data.summary, itemOpts{:});
            
            % ------ Selector for RTD distribution  -----------------------
            selector2D.dsel = tracerSelector(itemOpts{:});
            
            % create menu(s)
            % sub-menu for property display (3d and 2d)
            m1 = uiMenu('Title', 'Property display selection', 'Parent', d.Figure, ...
                        itemOpts{:}, 'items', {selector3D.psel, selector2D.ssel});
            % sub-menu for region selection
            m2 = uiMenu('Title', 'Region selection', 'Parent', d.Figure, ...
                        itemOpts{:}, 'items', {selector3D.wsel, selector3D.fsel});
            % sub-menu for diagnostics plots
            m3 = uiMenu('Title', 'Diagnostics', 'Parent', d.Figure, ...
                        itemOpts{:}, 'items', {selector2D.msel, selector2D.asel, selector2D.dsel});
            % main menu
            d.Menu = uiMenu('Title', 'Menu', 'Parent', d.Figure, ...
                itemOpts{:}, ...
                'items', {selector3D.tsel, m1, m2, m3});
            
            % set different bkcolor of each sub-menu ,select some light
            % colors
            tt = [251 180 76; 58 152 216; 42 187 155; 252 121 122]./255;
            tt = [tt; tatarizeMap];
            for k =1:numel(d.Menu.items)
                d.Menu.items{k}.BackgroundColor = tt(k,:);
                try
                    d.Menu.items{k}.titleColor = 'w';
                catch
                    d.Menu.items{k}.ForegroundColor = 'w';
                end
                for l=1:numel(d.Menu.items{k}.Children)
                    d.Menu.items{k}.Children(l).BackgroundColor = ...
                        .4*tt(k,:) + 0.6*[1 1 1];
                    for m=1:numel(d.Menu.items{k}.Children(l).Children)
                        d.Menu.items{k}.Children(l).Children(m).BackgroundColor = ...
                            .2*tt(k,:) + 0.8*[1 1 1];
                    end
                end
            end
            

            % axes
            d.Axes3D  = axes('Parent', d.Figure, 'Units', 'pixels', 'CLimMode', 'manual');
            d.Axes2DL = axes('Parent', d.Figure, 'Units', 'pixels');
            d.Axes2DR = axes('Parent', d.Figure, 'Units', 'pixels');
            [d.Axes2DL, d.Axes2DR] = addAxesContextMenu(d.Axes2DL, d.Axes2DR);
            
            
            % collapse submenues whose function do not yet make sense
            selector3D.wsel.collapse  = 1;
            for  f = reshape(fieldnames(selector2D),1, [])
               selector2D.(f{1}).collapse = 1;
            end
            
            % set figure callback functions:
            d.Figure.SizeChangedFcn        = @d.layout;
            d.Figure.WindowButtonDownFcn   = @d.resize;
            d.Figure.WindowButtonMotionFcn = @d.setMousePointer;
            d.Figure.WindowScrollWheelFcn  = @d.scrollMenu;
            
            % set zoom/pan/rotate
            d.interactiveModes           = setInteractiveModes(d.Axes3D);
            % ------ Show model with static property ----------------------
            d.Patch = cellDataPatch(d.G, d.Data.static(1).values, ...
               'Parent', d.Axes3D, 'EdgeColor', [.4 .4 .4], ...
               'EdgeAlpha', 1, 'BackFaceLighting', 'lit');
            d.Patch = addPatchContextMenu(d.Patch);
            
            d.Figure.CurrentAxes = d.Axes3D;
            d.outlineGrid = plotGrid(d.G, 'FaceColor', 'none', 'EdgeAlpha', 0.15, 'EdgeColor', [.4 .4 .4]);
            axis(d.Axes3D, 'tight', 'vis3d', 'off');
            d.Axes3D.ZDir = 'reverse';
            view(d.Axes3D, 3);
            daspect(d.Axes3D, [1 1 .2]);
            
            %return
            % ------ Add colorbar with histogram --------------------------
            d.colorBar = colorbar(d.Axes3D, 'WestOutside', 'AxisLocation', 'in');
            set(d.colorBar, 'Position',[.24 .8 .01 .19], 'Units', 'pixels');
            d.colorHAx = axes('Position',[.255 .8 .03 .19], 'Units', 'pixels', 'XDir', 'reverse');
            d.updateColorHist();
            
            % ------ Construct 2D plot axes -------------------------------
            axis(d.Axes2DL,'off');
            axis(d.Axes2DR,'off');
            
            % ------ Switch back to 3D plot axes and plot wells -----------
            % add extra toolbar stuff
            d.Figure.CurrentAxes = d.Axes3D;
            d.camlight    = light('Position', [1 1 -1], 'Style', 'infinite');
            d.camlight.Visible = 'off';
            d = addExtraTools(d);

            % finally wells
            d.WellPlot = wellPlotHandle(d.G, d.Data.states{1}.wellSol, ...
               'Visible', 'off', 'Parent', d.Axes3D);
            
            % ------ Set callbacks for 3D axes ----------------------------
            selector3D.tsel.Callback = @(src, event) d.tStepCallback(src, event, selector3D, selector2D);
            selector3D.psel.Callback = @(src, event) d.displayPropCallback(src, event, selector3D);
            selector3D.wsel.Callback = @(src, event) d.interactionRegionCallback(src, event, selector2D, selector3D);
            selector3D.fsel.Callback = @(src, event) d.filterPropCallback(src, event, selector3D);
            
            % ------ Set callbacks for 2D axes ----------------------------
            selector2D.msel.Callback      = @(src, event)d.measureCallback(src, event, selector2D, selector3D);
            selector2D.asel.Callback      = @(src, event)d.allocationCallback(src, event, selector2D, selector3D);
            %selector2D.dtsel.Callback     = @(src, event)tracerDistCallback(d, src, event, selector2D, selector3D);
            selector2D.ssel.Callback      = @(src, event)summaryCallback(d, src, event, selector2D);
            selector2D.ssel.regCallback   = @(src, event)selectWellsForSummary(d, src, event, selector2D, selector3D);
            selector2D.dsel.Callback      = @(src, event)d.distributionCallback(src, event, selector2D, selector3D);
           
            % ------ Do initial callback for 3D axes ----------------------
            selector3D.psel.Callback([], [])
             d.layout();
        end
            
        %% ---------- MAIN CALLBACKS --------------------------------------
        function tStepCallback(d, src, event, s3, s2)
            ts = s3.tsel.ix;
            d.currentDiagnostics = emptyDiagnostics(d);
            if isempty(ts) % disable all controls except
                s3.tsel.Enable = 'on';
                s3.wsel.Enable = 'off';   s3.wsel.collapse = 1;
                s2.msel.Enable = 'on';    s2.msel.collapse = 1;
                s2.asel.Enable = 'off';   s2.asel.collapse = 1;
                s2.dsel.Enable = 'off';   s2.dsel.collapse = 1;
                s2.dtsel.Enable = 'off'; s2.dtsel.collapse = 1;
                set([s3.psel.typePopup, s3.psel.propPopup], 'Enable', 'on');
                s3.psel.typeIx = 1;
                s3.psel.propIx = 4;
                s3.psel.propPopup.String = {d.Data.static.name};
                d.displayPropCallback(src, event, s3);
            else
                s3.wsel.Enable  = 'on';  s3.wsel.collapse = 0;
                s2.msel.Enable  = 'on';  s2.msel.collapse = 0;
                s2.dsel.Enable  = 'on'; 
                s2.asel.Enable  = 'on';
                s2.rsel.Enable  = 'on';
                s2.dtsel.Enable = 'on'; % Activate when this is implemented
                s3.wsel.collapse = 0;
                if numel(ts) == 1
                    % some stat-stuff
                end
                d.displayPropCallback(src, event, s3);
                d.interactionRegionCallback(src, event, s2, s3);
            end
        end
        % -----------------------------------------------------------------
        function displayPropCallback(d, src, event, s3) %#ok<*INUSL>
            switch s3.psel.typeIx
                case 1  % static property
                   displayVals = d.Data.static(s3.psel.propIx).values;
                   lims = d.Data.static(s3.psel.propIx).limits;
                case 2  % dynamic property
                    if ~isempty(s3.tsel.ix)
                        cval = d.Data.dynamic(s3.psel.propIx).values(:, s3.tsel.ix);
                        if numel(s3.tsel.ix) == 1
                            displayVals = cval;
                        else % multiple time steps
                            stat = s3.psel.statPopup.String{s3.psel.statIx};
                            displayVals = computeStatistic(cval, stat);
                        end
                        lims = d.Data.dynamic(s3.psel.propIx).limits;
                    else % no time-steps selected, can't disply dynamic property
                        s3.psel.typeIx = 1;
                        s3.psel.propPopup.String = {d.Data.static.name};
                        s3.psel.statPopup.String = {'n/a'};
                        s3.psel.statPopup.Enable = 'off';
                        return
                    end
                case 3   % diagnostics property
                    if ~isempty(s3.tsel.ix)
                        prop = s3.psel.propPopup.String{s3.psel.propIx};
                        [d, cval, lims, flag] = extractSelectedDiagnostics(d, prop, s3.tsel, s3.wsel);
                        if numel(s3.tsel.ix) == 1 || flag
                            displayVals = cval;
                            s3.psel.statIx = 1;
                            %s3.psel.statPopup.String = {'n/a'};
                            %s3.psel.statPopup.Enable = 'off';
                        else % multiple time steps
                            stat = s3.psel.statPopup.String{s3.psel.statIx};
                            displayVals = computeStatistic(cval, stat, prop); % include prop
                        end
                    else % no time-step selected, can't display diagnostic properties
                        s3.psel.typeIx = 1;
                        s3.psel.propPopup.String = {d.Data.static.name};
                        s3.psel.statPopup.String = {'n/a'};
                        s3.psel.statPopup.Enable = 'off';
                        return
                    end
                case 4   % computed property
                    if ~isempty(d.Props.computed.name)
                        displayVals = d.Data.computed(s3.psel.propIx).values;
                        lims        = d.Data.computed(s3.psel.propIx).limits;
                    else % no computed properties return to static
                        s3.psel.typeIx = 1;
                        s3.psel.propPopup.String = {d.Data.static.name};
                        return
                    end
            end
            isLog = false;
            if s3.psel.logSwitch
                [lims, displayVals, flag] = makeSafeForLog(lims, displayVals, 5);
                if ~flag % values not good for log-plot, reset switch
                    s3.psel.logSwitchBox.Value = 0;
                else
                    isLog = true;
                    [lims, displayVals] = deal(log10(lims), log10(displayVals));
                end
            end
            d.Patch.colorData = displayVals;
            d.Axes3D.CLim = lims;
            d.updateColorHist();
            d.updateColorBar(s3, lims, isLog);
        end        
        % -----------------------------------------------------------------        
        function filterPropCallback(d, src, event, s3)
            if ~strcmp(s3.fsel.playMode, 'stop')
                d.Patch.value = s3.fsel.maxValue;
            else
                if strncmp(src.Style,'check',5) && src.Value == 0 && strncmp(src.String,'Enable',6)
                    d.currentFilterRegion = [];
                    if ~isempty(d.currentInteractionRegion)
                        d.Patch.cells = d.currentInteractionRegion;
                    else
                        d.Patch.cells = 1:d.G.cells.num;
                    end
                    return
                end
                switch s3.fsel.typeIx
                    case 1  % static property
                        fval = d.Data.static(s3.fsel.propIx).values;
                    case 2  % dynamic property
                        if ~isempty(s3.tsel.ix)
                            fval = d.Data.dynamic(s3.fsel.propIx).values(:, s3.tsel.ix);
                            if numel(s3.tsel.ix) > 1 % multiple time steps
                                stat = s3.fsel.statPopup.String{s3.fsel.statIx};
                                fval = computeStatistic(fval, stat);
                            end
                        else % no time-steps selected, can't display dynamic property
                            s3.fsel.typeIx = 1;
                            s3.fsel.propPopup.String = {d.Data.static.name};
                            s3.fsel.statPopup.String = {'n/a'};
                            s3.fsel.statPopup.Enable = 'off';
                            updateSlideBars(s3.fsel, src, event)
                            fval = d.Data.static(s3.fsel.propIx).values;
                            d.Patch.colorData = fval;
                        end
                    case 3   % diagnostics property
                        if ~isempty(s3.tsel.ix)
                            prop = s3.fsel.propPopup.String{s3.fsel.propIx};
                            [d, fval] = extractSelectedDiagnostics(d, prop, s3.tsel, s3.wsel);
                            if numel(s3.tsel.ix) > 1 % multiple time steps
                                stat = s3.fsel.statPopup.String{s3.fsel.statIx};
                                fval = computeStatistic(fval, stat, prop); % include prop
                            end
                        else % no time-steps selected, can't display dynamic property
                            s3.fsel.typeIx = 1;
                            s3.fsel.propPopup.String = {d.Data.static.name};
                            s3.fsel.statPopup.String = {'n/a'};
                            s3.fsel.statPopup.Enable = 'off';
                            updateSlideBars(s3.fsel, src, event)
                            fval = d.Data.static(s3.fsel.propIx).values;
                            d.Patch.colorData = fval;
                        end
                    case 4 % computed prop
                        if ~isempty(d.Props.computed.name)
                            fval = d.Data.computed(s3.fsel.propIx).values;
                            %lims        = d.Data.computed(s3.psel.propIx).limits;
                        else % no computed properties return to static
                            s3.psel.typeIx = 1;
                            s3.psel.propPopup.String = {d.Data.static.name};
                            return
                        end
                end
                d.Patch.filterData = fval;
                d.currentFilterRegion = and(s3.fsel.minValue <= fval, s3.fsel.maxValue >= fval);
                if isempty(d.currentInteractionRegion)
                    d.Patch.cells = d.currentFilterRegion;
                else
                    d.Patch.cells =  and(d.currentFilterRegion, d.currentInteractionRegion);
                end
                d.updateColorHist();
            end
        end        
        % -----------------------------------------------------------------        
        function interactionRegionCallback(d, src, event, s2, s3)
            d.WellPlot.visibleInjectors = s3.wsel.injectorIx;
            d.WellPlot.visibleProducers = s3.wsel.producerIx;
            if isprop(src,'Style') && all(strncmp(src.Style,'checkbox',6)) 
               return;
            end
            % clear current diagnostics
            % this should only be done if well-selection has
            % changed, not if threshold is changed -> split into two
            % callbacks
            d.currentDiagnostics = emptyDiagnostics(d);
            [d, val] = extractSelectedDiagnostics(d, 'Tracer product', s3.tsel, s3.wsel);
            if numel(s3.tsel.ix) > 1 % do mean
                val = computeStatistic(val, 'mean');
            end
            d.currentInteractionRegion = val >= s3.wsel.threshold;
            if isempty(d.currentFilterRegion)
                d.Patch.cells = d.currentInteractionRegion;
            else
                if isprop(src,'Style') && all(strncmp(src.Style,'listbox',6))
                   d.filterPropCallback(src, event, s3);
                else
                   d.Patch.cells = ...
                      and(d.currentFilterRegion, d.currentInteractionRegion);
                end
            end
            if s3.psel.typeIx == 3 % reset display
                d.displayPropCallback(src, event, s3);
            end
            d.updateColorHist();
            % when interaction region has changed, set small axis color to
            % gray to indicate not up-to-date
            d.Axes2DL.Color = [.85 .85 .85];%d.Figure.Color;
            d.Axes2DR.Color = [.85 .85 .85];%d.Figure.Color;
            if s2.ssel.regionSwitch.Value == 1
               s2.ssel.regCallback(src, event)
            end
        end        
        % -----------------------------------------------------------------        
        function measureCallback(d, src, event, s2, s3)
           if s2.msel.panelNo==1
              ax = d.Axes2DL;
              d.resetSelectors(s2, 'msel','leftIx');
           else
              ax = d.Axes2DR;
              d.resetSelectors(s2, 'msel', 'rightIx');
           end
           cla(ax, 'reset');
           resetValue = false;
           switch src.Value
              case 1
                 axis(ax,'off');
              case 2
                 showWellCommunication(d, ax, s3.wsel.getCommunicationStrength());
              case 3
                 if isempty(s3.tsel.ix)
                    axis(ax,'off'); resetValue = true;
                 else
                    d.showFPhi(ax, s3.tsel, s3.wsel);
                 end
              case 4
                 if isempty(s3.tsel.ix)
                    axis(ax,'off'); resetValue = true;
                 else
                    d.showSweep(ax, s3.tsel, s3.wsel);
                 end
              case 5
                 if isempty(s3.tsel.ix)
                    axis(ax,'off'); resetValue = true;
                 else
                    d.showLorenz(ax, s3.tsel, s3.wsel);
                 end
              otherwise
                 disp('functionality not implemented yet');
           end
           if resetValue
              if s2.msel.panelNo==1
                 s2.msel.leftIx = 1;
              else
                 s2.msel.rightIx = 1;
              end
           end
        end        
        % -----------------------------------------------------------------        
        function allocationCallback(d, src, event, s2, s3)
            if s2.asel.panelNo==1
              ax = d.Axes2DL;
              d.resetSelectors(s2, 'asel', 'leftIx');
           else
              ax = d.Axes2DR;
              d.resetSelectors(s2, 'asel', 'rightIx');
           end
           cla(ax, 'reset');
           showAllocation(d, src, ax, s2, s3)
        end       
        % -----------------------------------------------------------------   
        function distributionCallback(d, src, event, s2, s3)
            [dsel, tsel, wsel] = deal(s2.dsel, s3.tsel, s3.wsel);
            if isempty(tsel.ix) || dsel.extendTime <= 0, return, end
            if s2.dsel.panelNo==1
                [ax, ix] = deal(d.Axes2DL, dsel.leftIx);
            else
                [ax, ix] = deal(d.Axes2DR, dsel.rightIx);
            end 
            cla(ax, 'reset');
            if numel(wsel.injectorIx) ~= 1 || numel(tsel.ix) ~= 1
                axis(ax, 'off')
                text(0,0,'Please select{\bf one} injector and{\bf one} time step.,','Parent' ,ax);
            else
                [iIx, pIx] = deal(s3.wsel.injectorIx, s3.wsel.producerIx);
                switch ix
                    case 1
                        axis(ax, 'off')
                    case 2 % estimate
                        %for k = 1:numel(pIx)
                            dist = estimateRTD(d.Gs.PORV, d.Data.diagnostics(tsel.ix).D, ...
                                                    d.Data.diagnostics(tsel.ix).WP, ...
                                                    'injectorIx', iIx, 'producerIx', pIx);
                            plot(ax, dist.t/year, dist.values, 'LineWidth', 2);
                        %end
                    case 3 % compute
                            dist = computeRTD(d.Data.states{tsel.ix}, d.Gs, d.Gs.PORV, ...
                                            d.Data.diagnostics(tsel.ix).D, d.Data.diagnostics(tsel.ix).WP, ...
                                            d.Data.states{tsel.ix}.wellSol, ...
                                            'injectorIx', iIx, 'producerIx', pIx);
                            plot(ax, dist.t/year, dist.values, 'LineWidth', 2);
                end
                if ix > 1
                    ylabel(ax, 'Tracer rate');
                    wn = arrayfun(@(x)x.label.String, d.WellPlot.producers(pIx), 'UniformOutput', false);
                    legend(ax, wn, 'Location','northeast', 'Interpreter', 'none')
                    set(ax, 'FontSize', 10)
                    set(ax, 'XLim', [0, s2.dsel.extendTime]);
                end
            end
        end
        % -----------------------------------------------------------------
        function summaryCallback(d, src, event, s2, ax)
            if nargin < 5
                if s2.ssel.panelNo == 1
                    ax = d.Axes2DL;
                else
                    ax = d.Axes2DR;
                end
            end
            [nms, prps] = deal(s2.ssel.curNames, s2.ssel.curProps);
            
            if isempty(prps)
                %axis(ax,'off');
            else
                if nargin < 5
                    cla(ax, 'reset');
                end
                leg = {};
                hold(ax, 'on')
                for k = 1:numel(nms)
                    for l = 1:numel(prps)
                        if any(strcmp(prps{l}, d.Data.summary.getKws(nms{k})))
                            plot(ax, s2.ssel.time, d.Data.summary.get(nms{k}, prps{l}, ':'), 'LineWidth', 2);
                            leg = [leg, {[nms{k},' - ', prps{l},' [', strtrim(d.Data.summary.getUnit(nms{k}, prps{l})), ']']}]; %#ok
                        end
                    end
                    if ~isempty(leg)
                        %xlabel('time [years]')
                        legend(ax, leg, 'Interpreter', 'none')
                        xtickangle(ax, 30)
                    end
                end
            end
        end
        % -----------------------------------------------------------------
        function selectWellsForSummary(d, src, event, s2, s3)
            prd = s3.wsel.prodSelector.String(s3.wsel.producerIx);
            inj = s3.wsel.injSelector.String(s3.wsel.injectorIx);
            s2.ssel.setWellSubset([prd(:);inj(:)]);
        end
        % -----------------------------------------------------------------
        function startPlayCallback(d, src, event)
            % set mouse pointer to wait-mode
            fsel = src.Parent.UserData;
            mtfnc  = d.Figure.WindowButtonMotionFcn;
            curpnt = d.Figure.Pointer; 
            d.Figure.WindowButtonMotionFcn = '';
            d.Figure.Pointer = 'watch'; 
            pause(0.001); % change pointer now 
            set(d.Axes3D, {'YLimMode', 'XLimMode', 'ZLimMode'}, {'manual', 'manual', 'manual'});
            d.Patch.logScale = fsel.logSwitch;
            d.Patch.playMode = true;
             % disable all exept playbar/slider
            cellfun(@(x)set(x, 'Enable', 'off'), d.Menu.items)
            set(fsel.playBar, 'Enable', 'on');
            fsel.maxSlider.Enable = 'on';
            %reset pointer
            d.Figure.Pointer = curpnt;
            d.Figure.WindowButtonMotionFcn = mtfnc;
        end
        % -----------------------------------------------------------------
        function stopPlayCallback(d, ~, ~)
            d.Patch.playMode = false;
            set(d.Axes3D, {'YLimMode', 'XLimMode', 'ZLimMode'}, {'auto', 'auto', 'auto'});
            % enable all
            cellfun(@(x)set(x, 'Enable', 'on'), d.Menu.items)
        end
        % -----------------------------------------------------------------
        %%  -------- LAYOUT CALLBACKS -------------------------------------
        function layout(d, ~, ~)
            [mw, ah, sp] = deal(d.layoutParams.menuWidth, d.layoutParams.distAxisHeight, d.layoutParams.itemSpace);
            fip    = d.Figure.Position;
            mPos   = d.Menu.Position;
            mPos   = [5, fip(4)-mPos(4)-1, mw, mPos(4)];
            mPos(3) = max(mPos(3), sp);
            aPos2D = [mw+sp, 2*sp, (fip(3)-mw-4*sp)/2, ah];
            aPos2D(4) = max(aPos2D(4), 0);
            aPos2DL = aPos2D; aPos2DL([1 3]) = aPos2DL([1 3]) + [sp -sp];
            aPos2DR = aPos2D; aPos2DR(1) = aPos2D(1)+aPos2D(3)+2*sp;
            aPos3D = [mw+sp, 2*sp+ah, fip(3)-mw-2*sp, fip(4)-3*sp-ah];
            cbh    = max(50, min(300, fip(4)-2*sp));
            cbw    = 27;
            %aPosCB = [mw+2*sp,       fip(4)-cbh-sp, cbw, cbh];
            aPosCB = [fip(3)-2*cbw-sp, fip(4)-cbh-sp, cbw, cbh];
            %aPosHA = [mw+2*sp+cbw+5, fip(4)-cbh-sp, cbw, cbh];
            aPosHA = [fip(3)-cbw-sp+5, fip(4)-cbh-sp, cbw, cbh];
            d.Menu.Position     = mPos;
            d.Axes2DL.Position  = aPos2DL;
            d.Axes2DR.Position  = aPos2DR;
            d.Axes3D.Position   = max(0,aPos3D);
            d.colorBar.Position = aPosCB;
            d.colorHAx.Position = aPosHA;
        end
        % -----------------------------------------------------------------
        function resize(d, src, ~)
            [mw, ah, sp] = deal(d.layoutParams.menuWidth, d.layoutParams.distAxisHeight, d.layoutParams.itemSpace);
            p = src.CurrentPoint;
            if p(1) <= mw % inside menu, deal with it there
                % do nothing
            elseif p(1) < mw+sp   % drag menu width
                motionFcn = d.Figure.WindowButtonMotionFcn;
                offset    = p(1)-mw;
                src.WindowButtonMotionFcn = {@dragHorisontally, offset};
                src.WindowButtonUpFcn     = {@clearMotionFcn, motionFcn};
            elseif (p(2)>ah+2*sp) && (p(2)<ah+3*sp)
                motionFcn = d.Figure.WindowButtonMotionFcn;
                offset    = p(2)-ah-sp;
                src.WindowButtonMotionFcn = {@dragVertically, offset};
                src.WindowButtonUpFcn     = {@clearMotionFcn, motionFcn};
            end
            
            function dragHorisontally(src, ~, offset)
                p = src.CurrentPoint - [offset, 0];
                if abs(mw - p(1)) > 5
                    d.layoutParams.menuWidth = max(p(1), 50);
                    d.Menu.Position(3) = d.layoutParams.menuWidth;
                    %d.layout();
                end
            end
            
            function dragVertically(src, ~, offset)
                p = src.CurrentPoint - [0 offset];
                if abs(ah - p(2)) > 2
                    d.layoutParams.distAxisHeight = max(p(2)-sp, sp);
                    d.layout();
                end
            end
            
            function clearMotionFcn(src, ~, motionFcn)
                src.WindowButtonMotionFcn = motionFcn;
                d.layout();
            end
        end
        %------------------------------------------------------------------
        function scrollMenu(d, src, event)
            mw = d.layoutParams.menuWidth;
            p = src.CurrentPoint;
            if p(1) < mw % over menu
                fp = d.Figure.Position;
                mp = d.Menu.Position; 
                if event.VerticalScrollCount > 0 
                    if mp(2) < 5
                        n = event.VerticalScrollAmount;
                        d.Menu.Position(2) = min(mp(2)+n*30, 5);
                    end
                elseif event.VerticalScrollCount < 0 
                    if mp(2)+ mp(4) > fp(4)-10
                        n = event.VerticalScrollAmount;
                        d.Menu.Position(2) = max(mp(2)-n*25, fp(4)-mp(4)-1);
                    end
                end
            end
        end        
        % -----------------------------------------------------------------
        function setMousePointer(d, ~, ~)
            [mw, ah, sp] = deal(d.layoutParams.menuWidth, d.layoutParams.distAxisHeight, d.layoutParams.itemSpace);
            p = d.Figure.CurrentPoint;
            if p(1) <= mw + sp
                if p(1) <= mw % inside menu panel
                    dy = arrayfun(@(x)x.Position(2), d.Menu.panel.Children)-p(2);
                    if any(and(dy > 0,dy <= 5))
                        set(d.Figure,'Pointer','top');
                    else
                        set(d.Figure,'Pointer','arrow');
                    end
                elseif ~any(structfun(@(x)strcmp(x.Enable, 'on'), d.interactiveModes)) % resize not allowed in interactive mode
                    set(d.Figure,'Pointer','left');
                else
                    set(d.Figure,'Pointer','arrow');
                end
            elseif (p(2)>ah+2*sp) && (p(2)<ah+3*sp) % between axes
                if ~any(structfun(@(x)strcmp(x.Enable, 'on'), d.interactiveModes))
                    set(d.Figure,'Pointer','top');
                else
                    set(d.Figure,'Pointer','arrow');
                end
            elseif ~any(structfun(@(x)strcmp(x.Enable, 'on'), d.interactiveModes))
                set(d.Figure,'Pointer','arrow');
            end
        end
        % -----------------------------------------------------------------
        %% --------- FUNCTIONS PLOTTING -----------------------------------
        function plotPie(d, ax, data, names, cmap)
           [data, pind] = sort(data, 'descend');
           dsum = sum(data);
           ix = find(cumsum(data) >= .95*dsum, 1, 'first');
           [data, pind] = deal(data(ix:-1:1), pind(ix:-1:1));
           h = pie(ax, max([data dsum-sum(data)]/dsum,eps));
           cmap = [cmap(pind,:); .8 .8 .8];
           for i=0:numel(pind)
              set(h(i*2+1),'FaceColor',.85*cmap(i+1,:)+[.15 .15 .15]);
           end
           legend(ax, [names(pind); 'others'], 'Location','EastOutside', 'Interpreter','none');
        end
        % -----------------------------------------------------------------
        function plotPLT(d, ax, wp, cmap)
           nseg = numel(wp.z);
           if nseg == 1  % Need to trick Matlab
              z = [wp.z; wp.z+1];
              a = [cumsum(wp.alloc,1); zeros(1,numel(wp.alloc))];
              bwidth = 2*numel(wp.z)+3;
              h=barh(ax, z, a,'stacked','BarWidth',bwidth, 'EdgeColor','none');
              args = {'YDir', 'YLim'};
           else
              z = flipud(wp.z);
              if any(diff(z)>=0), z=(nseg:-1:1).'; end
              a = cumsum(flipud(wp.alloc),1);
              if nseg<21
                 args = {'YDir', 'YLim'};
                 h=barh(ax, z, a,'stacked','BarWidth',1, 'EdgeColor','none');
              else
                 h=area(ax, z, a, eps,'EdgeColor','none'); axis tight;
                 view(90,-90);
                 args = {'XDir', 'XLim'};
              end
           end
           [zm, zM] = deal(min(z), max(z));
           for i=1:numel(h) % avoid using black bars
              set(h(i),'FaceColor',cmap(i,:))
           end
           set(ax,args{1},'reverse', args{2}, [zm zM + sqrt(eps)] + [-.1 .1]*(zM-zm));
        end
        % -----------------------------------------------------------------
        function plotPLT3D(d, ax, inj, reverseY, cmap)
           tmp = cellfun(@(x) x.alloc, inj, 'UniformOutput', false);
           alloc = cat(3, tmp{:});
           alloc = cumsum(alloc(end:-1:1,:,:), 1);
           z     = (size(alloc,1):-1:1).';
           hold(ax,'on');
           for i = 1:size(alloc,3)
              h=bar3h(ax, z, alloc(:,:,i),'stacked');
              for j=1:numel(h) 
                 set(h(j), 'xdata', get(h(j),'xdata')+i-1, ...
                    'FaceColor', cmap(j,:),'EdgeAlpha',.5);
              end
           end
           box(ax,'on'); grid(ax,'on'); view(ax,[140,20]); axis(ax,'tight');
           set(ax,'Xdir','reverse');
           if reverseY, set(ax,'Ydir','reverse'), end
        end
        % -----------------------------------------------------------------
        function showFPhi(d, ax, tsel, wsel)
           if isempty(tsel.ix), return, end
           m = getDynamicMeasures(d, tsel, wsel);
           
           nT = numel(tsel.ix);
           if nT < 8
              colM = get(gca,'ColorOrder');
           else
              colM = colorcube(nT+1);
           end
           for n=1:nT
              % Plot F-Phi for the whole regions
              col = colM(n,:);
              plot(ax,m.Phi,m.Ft(:,n),'LineWidth',2,'Color',col);
              if ~m.computePairs, continue; end
              
              % Plot F-Phi for individual wall pairs. For multiple time
              % steps, these lines are plotted using a lighter color. For a
              % single time step, we compute using different colors
              hold(ax,'on');
              np = size(m.F,2);
              if nT>1
                 col = repmat(.5*(col + [.6 .6 .6]), np, 1);
              elseif np<size(colM,1)
                 col = colM(2:end,:);
              else
                 col = colorcube(np+1);
              end
              for i=1:np
                 plot(ax,m.Phi,m.F(:,i,n),'Color',col(i,:), 'LineWidth',.5);
              end
           end
           hold(ax,'off');
           if ~isempty(m.wellName) 
              title(ax,['Well: ' m.wellName]);
              if nT<2
                 h=legend(ax,'region',m.names{:}); set(h,'FontSize',8)
              end
           end
        end
        % -----------------------------------------------------------------        
        function showSweep(d, ax, tsel, wsel)
           if isempty(tsel.ix), return, end
           m = getDynamicMeasures(d, tsel, wsel);
           
           nT = numel(tsel.ix);
           if nT < 8
              colM = get(gca,'ColorOrder');
           else
              colM = colorcube(nT+1);
           end
           for n=1:nT
              % Plot F-Phi for the whole regions
              col = colM(n,:);
              plot(ax,m.tDt(:,n),m.Ev,'LineWidth',2,'Color',col);
              if ~m.computePairs, continue; end
              
              % Plot F-Phi for individual wall pairs. For multiple time
              % steps, these lines are plotted using a lighter color. For a
              % single time step, we compute using different colors
              hold(ax,'on');
              np = size(m.F,2);
              if nT>1
                 col = repmat(.5*(col + [.6 .6 .6]), np, 1);
              elseif np<size(colM,1)
                 col = colM(2:end,:);
              else
                 col = colorcube(np+1);
              end
              for i=1:np
                 plot(ax,m.tD(:,i,n),m.Ev,'Color',col(i,:), 'LineWidth',.5);
              end
           end
           hold(ax,'off');
           if ~isempty(m.wellName) 
              title(ax,['Well: ' m.wellName]);
              if nT<2
                 h=legend(ax,'region',m.names{:}); set(h,'FontSize',8)
              end
           end
           hold(ax,'off'); axis(ax,'tight');
           set(ax,'XLim',[0 min(5,max(m.tDt(:)))]);
        end
        % -----------------------------------------------------------------  
        function showLorenz(d, ax, tsel, wsel)
           if isempty(tsel.ix), return, end
           
           m = getDynamicMeasures(d, tsel, wsel);
           
           axes(ax)
           nT = numel(tsel.ix);
           if ~m.computePairs
              h=bar(ax, m.LCt,'LineWidth',1, ...
                 'FaceColor',[.9 .9 .9],'EdgeColor',[0 .5 .8]);
              hold(ax,'on');
              text(h.XData, h.YData, cellstr(num2str(m.LCt(:),'%.4f')), ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment','bottom', 'FontSize', 8);
              
           elseif nT<4
              n = size(m.LC,1);
              h=bar(ax, m.LC,'LineWidth',1, ...
                 'FaceColor',[.9 .9 .9],'EdgeColor',[0 .5 .8]);
              grid(ax,'on'); hold(ax,'on');
              plot(ax,[.65 n+.35],[m.LCt m.LCt]','Color','r', 'LineWidth',2);
              if nT==1
                 text(h.XData, h.YData, cellstr(num2str(m.LC(:),'%.4f')), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment','bottom', 'FontSize', 8);
              else
                 yb = cat(1, h.YData);
                 xb = bsxfun(@plus, h(1).XData, [h.XOffset]');
                 text(xb(:),yb(:)-.05, cellstr(num2str(yb(:),'%.4f')), ...
                    'rotation', 90, 'horiz', 'right', 'FontSize',6);
              end
              set(ax,'XTickLabel',m.names,'FontSize',8);
              if n>8, set(ax,'XAxisLocation','top','XTickLabelRotation',45), end
           else
              n = size(m.LC,1);
              bar3(ax,m.LC); view(-80,20)
              hold(ax,'on')
              plot3(repmat(1:nT,2,1), repmat([.65 n+.35],nT,1)', ...
                 [m.LCt m.LCt]','Color','r','LineWidth',2);
              set(ax,'YTickLabel',m.names,'FontSize',8)
              if numel(m.names)>8, set(ax,'YTickLabelRotation',45), end
           end
           if ~isempty(m.wellName), title(['Well: ' m.wellName]); end
           hold(ax,'off'); axis(ax,'on', 'normal');
        end
        % -----------------------------------------------------------------
        %% --------- FUNCTIONS GENERAL ------------------------------------
        function resetSelectors(d, s, selector, field)
            flds = reshape(fieldnames(s),1, []);
            for f = flds(~strcmp(flds, selector))
                if isprop(s.(f{1}), field)
                    s.(f{1}).(field) = 1;
                end
            end
        end
        % -----------------------------------------------------------------            
        function updateColorHist(d)
           vals = d.Patch.colorData(d.Patch.cells);
           lims = d.Axes3D.CLim;
           [counts,bins] = hist(vals, linspace(lims(1), lims(2),25));
           visState = 'on';
           c = d.colorHAx.Children;
           if numel(c) == 1 && isprop(c, 'Visible')
               visState = d.colorHAx.Children.Visible;
           end
           barh(d.colorHAx, bins, counts,'hist');
           axis(d.colorHAx,'tight', 'off');
           h = findobj(d.colorHAx,'Type','patch');
           set(h,'FaceColor','none','EdgeColor',[.2 .2 .2], 'Visible', visState);
        end
        % -----------------------------------------------------------------
        function updateColorBar(d, s3, lims, islog)
            if nargin < 4
                islog = false;
            end
            cb = d.colorBar;
            % let first matlab choose ticks/labels
            if ~(s3.psel.typeIx == 3 && any(s3.psel.propIx == [7,8]))
                cb.TicksMode = 'auto';
                cb.TickLabelsMode = 'auto';
                cb.Limits = lims;
                if islog
                    p = cb.Ticks;
                    cb.TickLabels = cellfun(@(x)sprintf('%2.1e',10^x), num2cell(p), ...
                        'UniformOutput', false);
                end
                set(d.colorHAx.Children, 'Visible', 'on');
                colormap(d.Axes3D, 'default');
                %cb.Position(1) = d.layoutParams.menuWidth + 50;
            else
                cmap = tatarizeMap(max(size(d.Data.wellComunication,1), 64));
                if  s3.psel.propIx == 7 % sweep regions
                    ninj = numel(d.WellPlot.injectors);
                    cmap = cmap(2:ninj+1,:);
                    colormap(d.Axes3D, cmap);
                    cb.Ticks = (.5:ninj)/ninj;
                    cb.TickLabels = s3.wsel.injSelector.String;
                    cb.Limits = [0 1];
                elseif s3.psel.propIx == 8 % sweep regions
                    nprod = numel(d.WellPlot.producers);
                    cmap = cmap(2:nprod+1,:);
                    colormap(d.Axes3D, cmap);
                    cb.Ticks = (.5:nprod)/nprod;
                    cb.TickLabels = s3.wsel.prodSelector.String;
                    cb.Limits = [0 1];
                end
                set(d.colorHAx.Children, 'Visible', 'off');
                cb.TickLabelInterpreter = 'none';
                %cb.Position(1) = d.layoutParams.menuWidth + 75;
                % switch off log if selected
                if s3.psel.logSwitch == true
                    s3.psel.logSwitch = false;
                    s3.psel.logSwitchCallback(s3.psel.logSwitchBox, []);
                end
                %d.colorBar.Position = [360 500 30 350];
            end
            d.colorBar = cb;
        end
        % -----------------------------------------------------------------       
    end
end % --- END POSTPROCESSDIAGNOSICS ---------------------------------------

%% -- UTILITY FUNCTIONS ---------------------------------------------------
function s = emptyDiagnostics(d)
s = struct('name', '', 'values', []);
str = d.Props.diagnostics.name;
s = repmat(s, 1, numel(str));
for k = 1:numel(str)
    s(k).name = str{k};
end
end
% -----------------------------------------------------------------
function [lims, vals, flag] = makeSafeForLog(lims, vals, oom)
if lims(2) <= 0
    flag = false;
else
    flag = true;
    if lims(1) <= 0 || lims(1) < lims(2)*10^(-oom)
        lims(1) = lims(2)*10^(-oom);
        vals = max(vals, lims(1));
    end
end
end
% -----------------------------------------------------------------
function str = getTStepStrings(t_cur, t_prev)
formatOut = 'mmm dd, yyyy';
n = numel(t_cur);
str = cellfun(@(x)datestr(x, formatOut), mat2cell(t_cur, ones(1, n)), 'UniformOutput', false); 
if nargin > 1
    strp = cellfun(@(x)datestr(x, formatOut), mat2cell(t_prev, ones(1, n)), 'UniformOutput', false);
    str  = cellfun(@(x,y)[x,' - ',y], strp, str, 'UniformOutput', false);
end
end
% -----------------------------------------------------------------
% context menu helpers
function varargout = addAxesContextMenu(varargin)
nax = nargin;
varargout = cell(1,nax);
for k = 1:nax
    ax = varargin{k};
    assert(isa(ax, 'matlab.graphics.axis.Axes'));
    
    m = uicontextmenu('Parent', ax.Parent);
    uimenu(m, 'Label', 'Export to new figure', 'Callback', @(src, event)copyToFig(src, event, ax))

    ax.UIContextMenu = m;
    varargout{k} = ax;
end
end
% -----------------------------------------------------------------
function varargout = addPatchContextMenu(varargin)
nax = nargin;
varargout = cell(1,nax);
for k = 1:nax
    input = varargin{k};
    inputIsPatch = true;
    p = [];
    if isa(input, 'matlab.graphics.primitive.Patch')
        p = input;
    elseif isa(input, 'cellDataPatch')
        p = input.patchMain;
        inputIsPatch = false;
    end
    
    varargout{k} = input;
    
    if ~isempty(p)
        m = uicontextmenu('Parent', p.Parent.Parent);
        uimenu(m, 'Label', 'Export to new figure', 'Callback', @(src, event)copyToFig(src, event, p.Parent))
        if inputIsPatch
            varargout{k}.UIContextMenu = m;
        else
            varargout{k}.patchMain.UIContextMenu = m;
        end
    else
        warning('Unexpteted input: %s' ,class(input));
    end
end
end
% -----------------------------------------------------------------
function copyToFig(src, event, ax) %#ok
    f_new = figure;
    new = copyobj([ax, ax.Legend], f_new);
    ax  = findobj(new, 'Type', 'axes');
    drawnow
    set(ax, 'FontSize', 12, 'Units', 'normalized');
    set(ax, 'OuterPosition', [0 0 1 1])
end
% -----------------------------------------------------------------
% Toolbar customization
function d = addExtraTools(d)
tb  = findall(d.Figure,'Type','uitoolbar');
gh = uitoggletool(tb,'CData', getIcon('grid_outline.png'), 'Separator', 'on', 'State', 'on', ...
                   'TooltipString', 'Show/hide outline grid');
gh.OnCallback = @(src, event)outlineOnOff(src, event, d, 'on');
gh.OffCallback = @(src, event)outlineOnOff(src, event, d, 'off');

lh = uitoggletool(tb,'CData',getIcon('light.png'), 'Separator', 'off', 'State', 'off', ...
                   'TooltipString', 'Enable/disable light');
lh.OnCallback = @(src, event)lightOnOff(src, event, d, 'on');
lh.OffCallback = @(src, event)lightOnOff(src, event, d, 'off');  


aspFac = 1.2;
zp = uipushtool(tb,'CData',getIcon('z_plus.png'), 'Separator', 'on', ...
                'TooltipString','Increase vertical aspect ratio');
zp.ClickedCallback = @(src, event)changeVerticalAspectRatio(src, event, d, 1./aspFac);

zm = uipushtool(tb,'CData',getIcon('z_minus.png'), 'Separator', 'on', ...
                'TooltipString','Decrease vertical aspect ratio');
zm.ClickedCallback = @(src, event)changeVerticalAspectRatio(src, event, d, aspFac);

end
function outlineOnOff(~, ~, d, st)
d.outlineGrid.Visible = st;
end
function lightOnOff(~, ~, d, st)
if isvalid(d.camlight)
    d.camlight.Visible = st;
else
    d.camlight = camlight;
end
end
function changeVerticalAspectRatio(~, ~, d, fac)
    d.Axes3D.DataAspectRatio(3) = fac*d.Axes3D.DataAspectRatio(3);
end