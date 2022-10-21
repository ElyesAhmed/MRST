classdef WellboreModel < WrapperModel
    
    properties
        
        wells
        groups
        water, oil, gas
        thermal = false
        
    end
   
    methods
        
        %-----------------------------------------------------------------%
        function model = WellboreModel(rmodel, trajectories, varargin)
            
            opt = struct('names' , {{}}, ...
                         'WI'    , []  , ...
                         'groups', []);
            [opt, varargin] = merge_options(opt, varargin{:});
            % Get reservoir model subgrid consisting of all perforated cells
            if (iscell(trajectories) && isnumeric(trajectories{1})) ...
                    || isstruct(trajectories)
                % Well trajectories given as 
                [W, G] = processWellboreTrajectories( ...
                            rmodel.G, rmodel.rock, trajectories, ...
                            'names', opt.names, varargin{:});
            elseif iscell(trajectories) && isfield(trajectories{1}, 'cells')
                % Well trajectories given as grids
                error('Not supported yet')
            end
            
            % Get subset of rock
            rock = extractSubrock(rmodel.rock, G.cells.global);
            % Update reservoir model operators
            rmodel = rmodel.removeStateFunctionGroupings();
            rmodel.G = G; rmodel.rock = rock;
            rmodel = rmodel.setupOperators(G, rock);
            rmodel.useCNVConvergence = false;
            rmodel.nonlinearTolerance = 1e-3;
            phases = rmodel.getPhaseNames();
            for ph = phases
                rmodel.fluid.(['kr', ph']) = @(s) s;
            end
            
            % Construct WrapperModel
            model = model@WrapperModel(rmodel);
            model.wells = W;
            % Process groups
            if ~isempty(opt.groups)
                model.groups = model.processGroups(opt.groups);
            end

            % Set trajectories, default refDepths and names
            refDepth = zeros(model.numWells(), 1);
            names    = cell(model.numWells,1);
            for i = 1:model.numWells()
                cells       = G.cells.wellNo == i;
                refDepth(i) = min(G.cells.centroids(cells,3));
                names{i}    = sprintf('W%d', i);
            end
            model = merge_options(model, varargin{:});
            
            % Set phase and thermal flag
            model.water = model.parentModel.water;
            model.oil   = model.parentModel.oil;
            model.gas   = model.parentModel.gas;
            if isprop(model.parentModel, 'thermal')
                model.thermal = model.parentModel.thermal;
            end
            
            % Set up operators
            model = model.setupOperators(G, rock);
            if ~isempty(opt.WI)
                model.parentModel.operators.WI = opt.WI(model.G.cells.order);
            end
            model.operators = model.parentModel.operators;

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function groups = processGroups(model, groups)
            
            ng = numel(groups);
            wnames = {model.wells.name};
            memberIx = cell(model.numGroups, 1);
            for i = 1:ng
                group = groups(i);
                mix = find(ismember(wnames, group.members));
                memberIx{i} = reshape(mix, [], 1);
            end
            [groups.memberIx] = deal(memberIx{:});
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function G = computeWellGridGeometry(model, G)
            
            is_fractured = isfield(G.cells, 'hybrid');
            
            % Add inlet cells as copies of the top cell
            topCells = G.cells.topCell;
            [G, wcno] = copyCells(G, G.cells.topCell, 'nnc', repmat(topCells, 1, 2));
            % Set vertical coordinate equal to refDepth
            refDepth = arrayfun(@(W) W.refDepth, model.wells);
            G.cells.centroids(wcno,3) = refDepth;
            G.cells.global = [G.cells.global; G.cells.global(topCells)];
            if is_fractured
                G.cells.hybrid = [G.cells.hybrid; zeros(numel(wcno),1)]; 
            end
            
            gcno = [];
            if ~isempty(model.groups)
                % For groups, we use the top cell of the first member in
                % the group
                gw1 = arrayfun(@(group) group.memberIx(1), model.groups);
                gw  = arrayfun(@(group) group.memberIx, model.groups, 'UniformOutput', false);
                % Make nncs
                nnc = [rldecode(gw1, cellfun(@numel, gw), 1), vertcat(gw{:})];
                nnc = wcno(nnc); if size(nnc,2) == 1; nnc = nnc'; end
                % Copy cells
                [G, gcno] = copyCells(G, wcno(gw1), 'nnc', nnc);
                % Set vertical coord equal to refDepth of first group well
                G.cells.centroids(gcno,3) = refDepth(gw1);
                G.cells.global = [G.cells.global; G.cells.global(gw1)];
                if is_fractured
                    G.cells.hybrid = [G.cells.hybrid; zeros(numel(gcno),1)]; 
                end
            end
            
            % Add heat transmissibility fields for nncs
            [G.nnc.transHr, G.nnc.transHf] = deal(nan(size(G.nnc.trans)));
  
            % Set flag indicating if cell is perforation (i.e., connected
            % to a reservoir cell). Otherwise, it is an inlet(/outlet) cell
            G.cells.type = zeros(G.cells.num, 1);
            G.cells.type(wcno) = 1;
            G.cells.type(gcno) = 2;
            
            % Compute segment lengths and radii
            [lengthCell, lengthFace] = model.getSegmentLength(G);
            [radiusCell, radiusFace] = model.getSegmentRadius(G);
            
            % Get cell segment lengths and radii and set to G.cells for
            % easy access later
            G.cells.length = lengthCell;
            G.cells.radius = radiusCell;
            % Comute cross-sectional area of each cell segment
            area = pi.*G.cells.radius.^2;
            
            % Compute cell volumes
            G.cells.volumes = area.*G.cells.length;
            
            % Compute face area
            area = [nan; area];
            N    = G.faces.neighbors + 1;
            % Remember to rescale area-weighted normals
            G.faces.normals = G.faces.normals./G.faces.areas;
            % Compute as mean of well cross-section area of adjacent cells
            G.faces.areas   = mean(area(N),2,'omitnan');
            G.faces.normals = G.faces.normals.*G.faces.areas;
            
            G.connections = struct();
            G.connections.type = zeros(G.faces.num + numel(G.nnc.trans),1);
            G.connections.type((1:model.numWells) + G.faces.num) = 1;
            G.connections.type(G.faces.num+model.numWells+1:end) = 2;
            
            % Set face segment lengths and radii for easy access later
            G.faces.length = lengthFace;
            G.faces.radius = radiusFace;
            
        end
        
        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            
            % This will take care of parent model validation
            model = validateModel@WrapperModel(model, varargin{:});
            % Remove facility model
            model.parentModel.FacilityModel = [];
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            
            iic = model.getInletSegments();
            
            % Parent model state validation
            state = model.parentModel.validateState(state);
            % Add mass flux if not already set
            if ~isfield(state, 'massFlux')
                nf = nnz(model.parentModel.operators.internalConn);
                state.massFlux = zeros(nf, 1);
            end
            % Add bhp if not already set
            if ~isfield(state, 'bhp')
                bhp       = model.getProp(state, 'pressure');
                state.bhp = bhp(iic);
            end
            % Add surface rates if not already set
            rnames = model.getSurfaceRateNames();
            for rname = rnames
                if ~isfield(state, rname{1})
                    state.(rname{1}) = zeros(model.numWells() + model.numGroups(), 1);
                end
            end
            
            % Check if we have a thermal model
            if ~(isprop(model.parentModel, 'thermal') && model.parentModel.thermal)
                return
            end
            
            % Add bht if not already set
            if ~isfield(state, 'bht')
                bht       = model.getProp(state, 'temperature');
                state.bht = bht(iic);
            end
            % Add effect if not already set
            if ~isfield(state, 'effect')
                state.effect = zeros(model.numWells() + model.numGroups(), 1);
            end
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function rnames = getSurfaceRateNames(model)

            rnames = num2cell(model.parentModel.getPhaseNames());
            rnames = cellfun(@(name) ['q', name, 's'], rnames, 'UniformOutput', false); 

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)

            switch(lower(name))
                case {'massflux'}
                    fn = 'massFlux';
                    index = ':';
                case {'bhp'}
                    fn = 'bhp';
                    index = ':';
                case {'qws'}
                    fn = 'qWs';
                    index = ':';
                case {'qos'}
                    fn = 'qOs';
                    index = ':';
                case {'qgs'}
                    fn = 'qGs';
                    index = ':';
                case {'bht'}
                    fn = 'bht';
                    index = ':';
                case {'effect'}
                    fn = 'effect';
                    index = ':';
                otherwise
                    % Pass on to WrapperModel
                    [fn, index] = getVariableField@WrapperModel(model, name, varargin{:});
            end

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, varargin)

            % Parent model setup
            model = setupStateFunctionGroupings@WrapperModel(model, varargin{:});
            
            % Change flow discretization
            fd = model.parentModel.FlowDiscretization;
            % Replace component phase flux state function
            fd = fd.setStateFunction('ComponentPhaseFlux', WellboreComponentPhaseFlux(model));
            % Replace phase upwind flag state function
            fd = fd.setStateFunction('PhaseUpwindFlag', WellborePhaseUpwindFlag(model));
            % Add surface rate
            fd = fd.setStateFunction('SurfaceRate', WellboreSurfaceRate(model));
            % Add effect
            fd = fd.setStateFunction('Effect', WellboreEffect(model));
            % Set to model
            model.parentModel.FlowDiscretization = fd;
            
            % Change flow property functions
            fp = model.parentModel.FlowPropertyFunctions;
            % Add bottomhole pressure
            fp = fp.setStateFunction('BottomholePressure'   , WellboreBottomholePressure(model));
            % Add bottomhole temperature
            fp = fp.setStateFunction('BottomholeTemperature', WellboreBottomholeTemperature(model));
            % Set to model
            model.parentModel.FlowPropertyFunctions = fp;
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function model = setupOperators(model, G, rock, varargin)
            
            % Make wellbore grid and set to model
            nc0 = G.cells.num;
            GW = model.computeWellGridGeometry(G);
            
            % Make rock structure to represent wellbore
            rock = model.parentModel.rock;
            % Euivalent permeability from laminar pipe flow assumption.
            % This is now obsoloete due to dofs for the segment mass fluxes
            % + wellbore friction loss, but still kept here for reference.
            rock.perm = GW.cells.radius.^2/8;
            rock.poro = ones(GW.cells.num, 1);
            fnames = fieldnames(rock)';
            % Add values to rock fields for inlet cells
            cells = GW.cells.topCell;
            if ~isempty(model.groups)
                gw1 = arrayfun(@(group) group.memberIx(1), model.groups);
                cells = [cells; gw1];
            end
            for name = fnames
                v = rock.(name{1});
                if size(v,1) == nc0
                    rock.(name{1}) = [v; v(cells,:)];
                end
            end
            % Set rock and grid to model/parentModel
            model.parentModel.rock         = rock;
            [model.parentModel.G, model.G] = deal(GW);
            
            % Compute operators in standard way
            model.parentModel = model.parentModel.setupOperators(GW, rock, varargin{:});
            op = model.parentModel.operators;
            names = {'T', 'Thr', 'Thf'};
            [~, iif] = model.getInletSegments();
            for name = names
                if ~isfield(op, name{1}), continue; end
                fix = isnan(op.(name{1}));
%                 assert(~any(fix(~iif)), 'Computed non-finite transmissibilites');
                op.(name{1})(fix) = 0;
                op.([name{1}, '_all'])(op.internalConn) = op.(name{1});
            end
            model.parentModel.operators = op;
            
            % Set well indices
            WI = vertcat(model.wells.WI);
            WI = WI(model.G.cells.order);
            model.parentModel.operators.WI = WI;
            
            % Make mapping from interfaces to wells/groups
            nw = model.numWells(); ng = model.numGroups();
            mix = arrayfun(@(group) group.memberIx, model.groups, 'UniformOutput', false);
            ngw = cellfun(@numel, mix);
            ii = [(1:model.numWells)'; rldecode((1:numel(mix))', ngw, 1) + nw];
            jj = [(1:model.numWells)'; (1:sum(ngw))' + nw];
            FS = sparse(ii, jj, 1, nw + ng, nw + sum(ngw));
            
            model.parentModel.operators.fluxSum = @(v) FS*v;
            
            % Make mapping from wells to groups
            ii = [(1:model.numWells)'; rldecode((1:numel(mix))', ngw, 1) + nw];
            jj = [(1:model.numWells)'; vertcat(mix{:})];
            GS = sparse(ii, jj, 1, nw + ng, nw + sum(ngw));
            model.parentModel.operators.groupSum = @(v) GS*v;
            
            if model.thermal
                WIth = vertcat(model.wells.WIth);
                WIth = WIth(model.G.cells.order);
                model.parentModel.operators.WIth = WIth;
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [vars, names, origin] = getPrimaryVariables(model, state)
        % Get primary variables for model. These are the parent model
        % primary variables, in addition to interface mass fluxes
            
            % Get primary variable set of parent model
            [vars, names, origin] = model.parentModel.getPrimaryVariables(state);
            % Add mass flux
            vm = model.getProps(state, 'massFlux');
            vars = [vars, {vm}];
            names  = [names, 'massFlux'];
            origin = [origin, class(model)];

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function state = initStateAD(model, state, vars, names, origin)

            parent = strcmpi(origin, class(model.parentModel));
            state = model.parentModel.initStateAD(state, vars(parent), names(parent), origin(parent));

            vars   = vars(~parent);
            origin = origin(~parent); %#ok
            names  = names(~parent);
            
            for i = 1:numel(names)
                state = model.setProp(state, names{i}, vars{i});
            end
            state = model.initStateFunctionContainers(state);
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            
            % Get parent model equations
            [eqs, names, types, state] ...
                = model.parentModel.getModelEquations(state0, state, dt, drivingForces);
            % Set controls
            [eqs, state] = model.setControls(state, eqs, drivingForces);
            % Mass fluxes
            [fluxEqs, fluxNames, fluxTypes, state] = model.getMassFluxEquations(state);
            eqs   = [eqs  , fluxEqs  ];
            names = [names, fluxNames];
            types = [types, fluxTypes];
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getMassFluxEquations(model, state)
        % Equations realting segment mass fluxes to the pressure by means
        % of a wellbore friction model
        
            [~, roughness] = model.getSegmentRoughness();

            d = model.G.faces.radius.*2;
            l = model.G.faces.length;

            [v, rho, mu, pot, flag] = model.parentModel.getProps(state, ...
                                        'ComponentPhaseFlux'      , ...
                                        'Density'                 , ...
                                        'Viscosity'               , ...
                                        'PhasePotentialDifference', ...
                                        'PhaseUpwindFlag'         );
            rho{1} = model.parentModel.operators.faceAvg(rho{1});

            mu{1} = model.parentModel.operators.faceUpstr(flag{1}, mu{1});
            dp = wellBoreFriction(v{1}, rho{1}, mu{1}, d, l, roughness, 'massRate');

            bhp = model.getProp(state, 'bhp');
            eqs   = (pot{1} - dp)./max(max(value(bhp)), 1);
            eqs   = {eqs};
            names = {'flux'};
            types = {'face'};

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [eqs, state] = setControls(model, state, eqs, drivingForces)
        % Equations imposing well control for each well
        
            % Get driving forces and check that it's non-empty
            [wCtrl, gCtrl] = deal([]);
            hasWellCtrl  = isfield(drivingForces, 'W');
            hasGroupCtrl = isfield(drivingForces, 'groups');
            if hasWellCtrl , wCtrl = drivingForces.W;      end
            if hasGroupCtrl, gCtrl = drivingForces.groups; end
            if isempty(wCtrl) && isempty(gCtrl), return; end
            
            nw = model.numWells();
            ng = model.numGroups();
            [iic, iif] = model.getInletSegments(); ic = find(iic);
            
            ctrl0 = struct('type', 'none', 'val', nan, 'T', nan);
            if ~hasWellCtrl , wCtrl = repmat(ctrl0, nw, 1) ; end
            if ~hasGroupCtrl, gCtrl = repmat(ctrl0, ng, 1); end
        
            ctrlTarget = [vertcat(wCtrl.val); vertcat(gCtrl.val)];
            ctrlType   = vertcat({wCtrl.type}', {gCtrl.type}');
            
            % No control - well is controlled by group or group is not
            % actively controlling wells
            is_none  = strcmpi(ctrlType, 'none');
            % Well is controlled by its group
            is_group = strcmpi(ctrlType, 'group');
            % Treat wells controlled by a group as "no control"
            is_none = is_none | is_group;
            % Bottom-hole pressure control
            is_bhp  = strcmpi(ctrlType, 'bhp');
            % Surface total rates
            is_rate = strcmp(ctrlType, 'rate') | strcmpi(ctrlType, 'vrat');
            % Surface oil rates
            is_orat = strcmp(ctrlType, 'orat');
            % Surface water rates
            is_wrat = strcmp(ctrlType, 'wrat');
            % Surface gas rates
            is_grat = strcmp(ctrlType, 'grat');
            % Surface liquid rates (water + oil)
            is_lrat = strcmp(ctrlType, 'lrat');
            % Reservoir rates (at averaged conditions from previous step)
            is_resv = strcmp(ctrlType, 'resv');
            % Reservoir rates (at current conditions for each perf.)
            is_volume = strcmp(ctrlType, 'volume');
            % Check for unsupported well controls
            assert(~any(is_resv | is_volume), 'Not implemented yet');
            
            bhp = model.getProps(state, 'BottomholePressure');
            % Set BHP control equations
            eqs{1}(ic(is_bhp)) = ctrlTarget(is_bhp) - bhp(is_bhp);

            % Set rate control equations
            nph    = model.parentModel.getNumberOfPhases();
            phases = model.parentModel.getPhaseNames();
            is_surface_rate = false(nw + ng, 1);
            % Find rate targets
            for ph = 1:nph
                switch phases(ph)
                    case 'W'
                        act = is_rate | is_wrat | is_lrat;
                    case 'O'
                        act = is_rate | is_orat | is_lrat;
                    case 'G'
                        act = is_rate | is_grat;
                end
                is_surface_rate(act) = true;
            end
            
            % Convert surface volume rates to mass rates
            rhoS = model.getProps(state, 'SurfaceDensity');
            rhoS = rhoS{1}(iic);
            Q = ctrlTarget(is_surface_rate).*rhoS(is_surface_rate);
            % Add in as source in corresponding well cells
            eqs{1}(ic(is_surface_rate)) = eqs{1}(ic(is_surface_rate)) - Q;
   
            if model.thermal
                
                % Get properties
                [bht, p, T, hr, qh, qr, rho] = model.getProps(state, ...
                    'BottomholeTemperature', 'pressure', 'T', ...
                    'enthalpy', 'HeatFlux', 'massFlux', 'Density');
                % Compute total well rates
                qt = model.parentModel.operators.fluxSum(qr(iif));
                qt = sum(value(qt), 2);
                % Identify injection wells/groups
                is_inj  = qt > 0 | (is_surface_rate & ctrlTarget > 0);
                is_inj  = is_inj & ~is_none;
                % Compute temperature in production well cells based on
                % conservation of energy
                qh = model.parentModel.operators.Div(qh);
                qr = model.parentModel.operators.Div(qr);
                is_zero = abs(qr) < 10*eps;
                h = abs(qh./qr);
                h(is_zero) = hr(is_zero);
                Cp = model.parentModel.fluid.CpW(p, T);
                T = (h - p./rho{1})./Cp;
                T = T(iic);
                % Set well temperature targets
                ctrlT = [vertcat(wCtrl.T); vertcat(gCtrl.T)];
                assert(~any(isnan(value(ctrlT(is_inj)))), ...
                    'Temperature in active well/group not set!');
                ctrlT(isnan(ctrlT)) = mean(ctrlT(~isnan(ctrlT)));
                Tw = ctrlT.*is_inj + T.*~is_inj;
                % Set equations
                cells = ~is_none & ~is_zero(iic);
                eqs{2}(ic(cells)) = (Tw(cells) - bht(cells));
                
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, forces)

            parent = strcmpi(problem.types, 'cell');

            p = problem; p.primaryVariables = p.primaryVariables(parent);
            [state, report] = updateState@WrapperModel(model, state, p, dx, forces);

            for i = (nnz(parent) + 1):numel(problem.primaryVariables)
                 p = problem.primaryVariables{i};
                 % Update the state
                 state = model.updateStateFromIncrement(state, dx, problem, p);
            end

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%        
        function [model, state] = prepareReportstep(model, state, state0, dt, drivingForces) %#ok
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            
            [convergence, values, names] = checkConvergence@WrapperModel(model, problem, varargin{:});
            
%             dt = problem.dt;
%             Q  = abs(value(problem.state.inletFlux{1}));
% %             Q = max(Q);
%             Q = Q(model.G.cells.wellNo);
%             Qh = abs(value(problem.state.inletHeatFlux));
% %             Qh = max(Qh);
%             Qh = Qh(model.G.cells.wellNo);
%             [mass, energy] = model.parentModel.getProps(problem.state, 'ComponentTotalMass', 'TotalThermalEnergy');
%             
%             pad = zeros(model.numWells() + model.numGroups() ,1);
%             
%             scaleMass = max(value(mass{1})./dt, [Q; pad]);
%             scaleEnergy = max(value(energy)./dt, [Qh; pad]);
%             
%             values(1) = norm(value(problem.equations{1}./scaleMass), inf);
%             values(2) = norm(value(problem.equations{2}./scaleEnergy), inf);
%             convergence(1:2) = values(1:2) < model.parentModel.nonlinearTolerance;
%             convergence(1:2) = values(1:2) < 1e-2;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            
            [state, report] = model.parentModel.updateAfterConvergence(state0, state, dt, drivingForces);
            
            [state.bhp, state.bht, qr] = model.getProps(state, ...
                'BottomholePressure', 'BottomholeTemperature', 'massFlux');

            [qh, rhoS, flag] = model.parentModel.getProps(state, ...
                'HeatFlux', 'SurfaceDensity', 'PhaseUpwindFlag');

            [iic, iif] = model.getInletSegments();
            % Upstream-wheight density
            rhoSf = cellfun(@(flag, rho) ...
                model.parentModel.operators.faceUpstr(flag, rho), ...
                flag, rhoS, 'UniformOutput', false);
            qs = model.parentModel.operators.groupSum(qr(iif)./rhoSf{1}(iif));
            qh = model.parentModel.operators.groupSum(qh(iif));
            % Set surface rate and effect
            state.qWs    = qs;
            state.effect = qh;
            
            [wCtrl, gCtrl] = deal([]);
            hasWellCtrl  = isfield(drivingForces, 'W');
            hasGroupCtrl = isfield(drivingForces, 'groups');
            if hasWellCtrl , wCtrl = drivingForces.W;      end
            if hasGroupCtrl, gCtrl = drivingForces.groups; end
            
            nw = model.numWells();
            ng = model.numGroups();
            
            ctrl0 = struct('type', 'none', 'val', nan, 'T', nan);
            if ~hasWellCtrl , wCtrl = repmat(ctrl0, nw, 1); end
            if ~hasGroupCtrl, gCtrl = repmat(ctrl0, ng, 1); end
        
            ctrlType = vertcat({wCtrl.type}', {gCtrl.type}');
            is_none  = strcmpi(ctrlType, 'none');
            
            CpW  = model.parentModel.fluid.CpW(state.bhp, state.bht);
            rhoW = model.parentModel.fluid.rhoW(state.bhp, state.bht);
            h = state.effect./(state.qWs.*rhoS{1}(iic));
            T = (h - state.bhp./rhoW)./CpW;
            % Set bottom-hole temperature in inactive groups
            state.bht(is_none) = T(is_none);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function n = numWells(model)
            
            n = numel(model.wells);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function n = numGroups(model)
            
            n = numel(model.groups);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function cells = getGlobalWellCells(model)
        % Get cell numbering of perforated cells in the reservoir grid
        
            cells = model.G.cells.global(model.G.cells.type == 0);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [iiCell, iiFace] = getInletSegments(model, G)
            
            if nargin < 2, G = model.G; end
            iiCell = G.cells.type ~= 0;
            iiFace = any(iiCell(model.parentModel.operators.N),2);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [rc, rf] = getSegmentRadius(model, G)
            
            if nargin < 2, G = model.G; end
            
            % Get radius in each cell
            rc = vertcat(model.wells.r);
            rc = rc(G.cells.order);
            
            % Add radii for inlet well cells
            rc = [rc; rc(G.cells.topCell)];
            
            if ~isempty(model.groups)
            % Add radii for inlet group cells
                gw1 = arrayfun(@(group) group.memberIx(1), model.groups);
                rc = [rc; rc(G.cells.topCell(gw1))];
            end
            
            % Face segment radii computed as mean of adjacent cells
            N = G.faces.neighbors; N = N(all(N>0,2), :);
            N = [N; G.nnc.cells];
            rf = mean(rc(N),2);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [lc, lf] = getSegmentLength(model, G)
            
            if nargin < 2, G = model.G; end
            
            % Compute lenght of segments within each cell
            lc = arrayfun(@(W) sqrt(sum(W.trajectory.vec.^2, 2)), ...
                        model.wells, 'UniformOutput', false);
            lc = cell2mat(lc);
            lc = lc(G.cells.order);

            % Compute lenght of segments for each inlet cell
            cells = G.cells.topCell;
            if model.numGroups() > 0
                % Get cells that group cells were copied from
                gw1 = arrayfun(@(group) group.memberIx(1), model.groups);
                cells = [cells; gw1];
            end
            lc = [lc; lc(cells)];
            
            % Face segment lengths computed as mean of adjacent cells            
            N = G.faces.neighbors; N = N(all(N>0,2), :);
            N = [N; G.nnc.cells];
            lf = sum(lc(N).*0.5, 2);
%             if model.numGroups > 0
%                 x  = G.cells.centroids;
%                 li = sqrt(sum((x(G.nnc.cells(:,1), :) ...
%                                 - x(G.nnc.cells(:,2), :)).^2, 2)); 
%                 lf = [lf; li];
%             end

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [rc, rf] = getSegmentRoughness(model, G)
            
            if nargin < 2, G = model.G; end
             
            % Get wellbore rougness (currently assumes one value per well)
            r0 = arrayfun(@(W) sum(W.trajectory.roughness, 2), ...
                    model.wells, 'UniformOutput', false);
            rc = cell2mat(r0);
            assert(numel(rc) == model.numWells, ['WellboreModel ', ...
                'currently only supports one roughness value per well']);
            rc = rc(G.cells.wellNo);
            cells = G.cells.topCell;
            if model.numGroups() > 0
                % Get cells that group cells were copied from
                gw1 = arrayfun(@(group) group.memberIx(1), model.groups);
                cells = [cells; gw1];
            end
            rc = [rc; rc(cells)];
            
            % Face segment roughness computed as mean of adjacent cells
            N = G.faces.neighbors; N = N(all(N>0,2), :);
            N = [N; G.nnc.cells];
            rf = sum(rc(N).*0.5, 2);

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function names = getComponentNames(model)
        % Get names of components in parent model
            
            names = model.parentModel.getComponentNames();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function ncomp = getNumberOfComponents(model)
        % Get number of components in parent model
            
            ncomp = model.parentModel.getNumberOfComponents();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function nph = getNumberOfPhases(model)
        % Get number of phases in parent model
            
            nph = model.parentModel.getNumberOfPhases();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function active = getActivePhases(model)
        % Get active phases in parent model
            
            active = model.parentModel.getActivePhases();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function rhoS = getSurfaceDensities(model, varargin)
            
            rhoS = model.parentModel.getSurfaceDensities(varargin{:});
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function plotWellNetwork(model)
            nnc = model.G.nnc.cells(model.numWells + 1:end, :);
            nnc = max(nnc(:)) - nnc + 1;
            names = [{model.wells.name}, {model.groups.name}];
            names = fliplr(names);
            network = graph(nnc(:,1), nnc(:,2), [], names);
            
            plot(network, 'layout', 'layered', 'LineWidth', 2, 'MarkerSize', 8);
            
%             cc = network.conncomp;
%             hold on;
%             for i = 1:max(cc)
%                 wix = cc == cc(i);
%                 plot(network.subgraph(wix), 'layout', 'layered', 'LineWidth', 2);
%             end
        end
        %-----------------------------------------------------------------%

    end
    
end

function rock = extractSubrock(rock, cells)

    n = numel(rock.poro);
    names = fieldnames(rock);
    for name = names'
        v = rock.(name{1});
        if size(v,1) == n
            rock.(name{1}) = v(cells,:);
        end
    end

end