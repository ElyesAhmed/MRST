classdef TransportModelDG < TransportModel
    
    properties
        discretization
        limiters
        storeUnlimited
        dgVariables
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = TransportModelDG(parent, varargin)
           
            model = model@TransportModel(parent);
            model.discretization = [];
            model.dgVariables = {'s'};
            
            % Default limiters
            names    = model.dgVariables;
            limits   = {[0,1]};
            tol      = 0;
            limiters = [];
            limiters = addLimiter(limiters, ...
                                  'type'     , 'tvb' , ...
                                  'variables', names , ...
                                  'limits'   , limits, ...
                                  'tol'      , tol   );
            limiters = addLimiter(limiters, ...
                                  'type'     , 'scale', ...
                                  'variables', names  , ...
                                  'limits'   , limits , ...
                                  'tol'      , tol    );
            model.limiters       = limiters;
            model.storeUnlimited = false;
            
            [model, discretizationArgs] = merge_options(model, varargin{:});
            % Construct discretization
            if isempty(model.discretization)
                model.discretization = DGDiscretization(model.G, discretizationArgs{:});
            end
            % Assign discretization to parentModel
            model.parentModel.discretization = model.discretization;
            if strcmpi(model.formulation, 'totalSaturation') && ...
                    any(strcmpi(model.dgVariables, 's'))
                model.dgVariables{end+1} = 'sT';
            end
            % Get limiters
            for l = 1:numel(model.limiters)
                limiter = model.limiters(l);
                model.limiters(l).function = getLimiter(model, limiter.type);
            end
            % Set up DG operators
            model.parentModel.operators = setupOperatorsDG(model.discretization  , ...
                                                           model.parentModel.G   , ...
                                                           model.parentModel.rock);
            % Pressure is not solved with DG, make sure to don't store
            % things to state that should be recomputed in pressure step
            model.parentModel.outputFluxes         = false;
            model.parentModel.OutputStateFunctions = {};
        end
        
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            % Get variable fiels, check if it is dof
            isDof = numel(name) > 3 && strcmpi(name(end-2:end), 'dof');
            if isDof
                lookup = name(1:end-3);
            else
                lookup = name;
            end
            [fn, index] = getVariableField@TransportModel(model, lookup, varargin{:});
            if isDof && ~isempty(fn)
                fn = [fn, 'dof'];
            end
        end
        
        %-----------------------------------------------------------------%
        function groupings = getStateFunctionGroupings(model)
            groupings = model.parentModel.getStateFunctionGroupings();
        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            % Set degree in each cell
            state.degree = repmat(model.discretization.degree, model.G.cells.num, 1);
            % Well are treated as dG(0)
            wm = model.parentModel.FacilityModel.WellModels;
            for i = 1:numel(wm)
                state.degree(wm{i}.W.cells,:) = 0;
            end
            % Let parent model do its thing
            state = validateState@TransportModel(model, state);
            % Assign dofs
            state = assignDofFromState(model.discretization, state, model.dgVariables);
        end
        
        %-----------------------------------------------------------------%
        function [dofvars, dofnames, names, origin] = getPrimaryVariables(model, state)
            % Get primary variables
            [vars, names, origin] = model.parentModel.getPrimaryVariables(state);
            isParent = strcmp(origin, class(model.parentModel));
            vars = vars(isParent);
            names = names(isParent);
            dofnames = cellfun(@(bn) [bn, 'dof'], names, 'UniformOutput', false);
            origin = origin(isParent);
            dofvars = cell(1, numel(vars));
            isBO = strcmpi(origin, 'GenericBlackOilModel');
            for i = 1:numel(dofnames)
                [fn, ~] = model.getVariableField(dofnames{i}, false);
                if ~isempty(fn)
                    dofvars{i} = model.getProp(state, dofnames{i});
                elseif strcmpi(names{i}, 'x') && isBO(i)
                    if model.parentModel.water
                        sW = model.getProp(state, 'sW');
                    else
                        sW = deal(0);
                    end
                    [sG, sGdof] = model.getProps(state, 'sG', 'sGdof');
                    st = model.parentModel.getCellStatusVO(state,  1-sW-sG, sW, sG);
                    for j = 1:numel(st)
                        if numel(st{j}) == model.G.cells.num
                            st{j} = rldecode(st{j}, state.nDof, 1);
                        end
                    end
                    [rsdof, rvdof] = model.getProps(state, 'rsdof', 'rvdof');
                    xdof = st{1}.*rsdof + st{2}.*rvdof + st{3}.*sGdof;
                    dofvars{i} = xdof;
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function [state, names, origin] = getStateAD(model, state, init)
            if nargin < 3
                init = true;
            end
            
            if 0
                names = fieldnames(state);
                cells = rldecode((1:model.G.cells.num)', state.nDof, 1);
                for k = 1:numel(names)
                    name = names{k};
                    if numel(name) > 3 && strcmp(name(end-2:end), 'dof')
                        v   = model.getProp(state, name(1:end-3));
                        dof = model.getProp(state, name);
                        vm  = model.discretization.getCellMean(state, dof);
                        frac = v./vm;
                        frac(~isfinite(frac)) = 1;
                        dof = dof.*frac(cells,:);
                        state = model.setProp(state, name, dof);
                    end
                end
            end
            
            parent = model.parentModel;
            % Get the AD state for this model
            [basevars, basedofnames, basenames, baseorigin] = model.getPrimaryVariables(state);
            isParent = strcmp(baseorigin, class(parent));
            basevars = basevars(isParent);
            basedofnames = basedofnames(isParent);
            basenames = basenames(isParent);
            baseorigin = baseorigin(isParent);
            % Find saturations
            isS = false(size(basevars));
            nph = parent.getNumberOfPhases();
            phase_variable_index = zeros(nph, 1);
            for i = 1:numel(basevars)
                [f, ix] = model.getVariableField(basedofnames{i}, false);
                if strcmp(f, 'sdof')% || strcmpi(basedofnames{i}, 'xdof')
                    isS(i) = true;
                    phase_variable_index(ix) = i;
                end
            end
            
            % Figure out saturation logic
            isP    = strcmp(basedofnames, 'pressuredof');
            vars   = basevars;
            names  = basedofnames;
            origin = baseorigin;
            useTotalSaturation = strcmpi(model.formulation, 'totalSaturation') ...
                                    && sum(isS) == nph - 1;
            useTotalSaturation = useTotalSaturation ...
                || strcmpi(class(parent), 'GenericOverallCompositionModel');
            if useTotalSaturation
                % Replace pressure with total saturation
                replacement = 'sTdof';
                sTdof = model.getProp(state, replacement);
                % Replacing
                vars{isP} = sTdof;
                names{isP} = replacement;
                origin{isP} = class(model);
            else
                % Remove pressure and skip saturation closure
                vars = vars(~isP);
                names = names(~isP);
                origin = origin(~isP);
            end
            if init
                [vars{:}] = model.AutoDiffBackend.initVariablesAD(vars{:});
            end
            if useTotalSaturation
                basevars(~isP) = vars(~isP);
            else
                basevars(~isP) = vars;
            end
            
            [~ , xc, cNo     ] = model.discretization.getCubature(Inf, 'cell');
            [~ , xf, ~  , fNo] = model.discretization.getCubature(Inf, 'face');
            xf   = repmat(xf, 2, 1);
            N    = model.discretization.G.faces.neighbors;
            fcNo = [N(fNo,1); N(fNo,2)];
            xc = model.discretization.transformCoords(xc, cNo );
            xf = model.discretization.transformCoords(xf, fcNo);
            [psi_c, psi_f] = deal(model.discretization.basis.psi');
            for dofNo = 1:model.discretization.basis.nDof
                psi_c{dofNo} = psi_c{dofNo}(xc);
                psi_f{dofNo} = psi_f{dofNo}(xf);
            end
            state.psi_c = psi_c;
            state.psi_f = psi_f;
            
            [cellMean, cellVars, faceVars] = deal(cell(size(vars)));
            for i = 1:numel(vars)
                cellMean{i} = model.discretization.getCellMean(state, basevars{i});
                cellVars{i} = model.discretization.evaluateProp(state, basevars{i}, 'cell');
                faceVars{i} = model.discretization.evaluateProp(state, basevars{i}, 'face');
            end
            
            state = parent.initStateAD(state, cellMean, basenames, baseorigin);
            state.cells = (1:model.G.cells.num)';
            state = model.evaluateBaseVariables(state);
            
            state.wellStateDG = parent.initStateAD(state.wellStateDG, cellMean, basenames, baseorigin);
            
            parent.G.cells.num = numel(value(cellVars{1}));
            state.cellStateDG = parent.initStateAD(state.cellStateDG, cellVars, basenames, baseorigin);
            
            parent.G.cells.num = numel(value(faceVars{1}));
            state.faceStateDG = parent.initStateAD(state.faceStateDG, faceVars, basenames, baseorigin);

            if useTotalSaturation
                % Set total saturation as well
                sTdof       = vars{isP};
                state.sTdof = sTdof;
                % Evaluate at cell cubature points
                cellValue         = model.discretization.evaluateProp(state, sTdof, 'cell');
                state.cellStateDG = model.setProp(state.cellStateDG, 'sT', cellValue);
                % Evaluate mean
                cellMean          = model.discretization.getCellMean(state, sTdof);
                state.wellStateDG = model.setProp(state.wellStateDG, 'sT', cellMean);
                % Evaluate at face cubature points
                faceValue         = model.discretization.evaluateProp(state, sTdof, 'face');
                state.faceStateDG = model.setProp(state.faceStateDG, 'sT', faceValue);
                % Set mean in state
                state = model.setProp(state, 'sT', cellMean);
            end
        end
        
        function state = evaluateBaseVariables(model, state)
             
            [cellStateDG, faceStateDG, wellStateDG] = deal(state);
            % Evaluate basis functions at cubature points
            [~ , xc, cNo     ] = model.discretization.getCubature(Inf, 'cell');
            [~ , xf, ~  , fNo] = model.discretization.getCubature(Inf, 'face');
            xf   = repmat(xf, 2, 1);
            N    = model.discretization.G.faces.neighbors;
            fcNo = [N(fNo,1); N(fNo,2)];
            xc = model.discretization.transformCoords(xc, cNo );
            xf = model.discretization.transformCoords(xf, fcNo);
            [psi_c, psi_f] = deal(model.discretization.basis.psi');
            for dofNo = 1:model.discretization.basis.nDof
                psi_c{dofNo} = psi_c{dofNo}(xc);
                psi_f{dofNo} = psi_f{dofNo}(xf);
            end
            state.psi_c = psi_c;
            state.psi_f = psi_f;
            % Set cells/faces
            [~, ~, cells] = model.discretization.getCubature((1:model.G.cells.num)', 'cell');
            [~, ~, ~, faces] = model.discretization.getCubature(find(model.parentModel.operators.internalConn), 'face');
            fcells = [model.G.faces.neighbors(faces,1); model.G.faces.neighbors(faces,2)];
            cellStateDG.type   = 'cell';
            cellStateDG.cells  = cells;
            cellStateDG.fcells = fcells;
            cellStateDG.faces  = faces;
            wellStateDG.type   = 'cell';
            cellStateDG.cells  = cells;
            faceStateDG.type   = 'face';
            faceStateDG.cells  = fcells;
            faceStateDG.faces  = faces;
            
            names = fieldnames(state);
            for k = 1:numel(names)
                name = names{k};
                [cellValue, cellMean, faceValue] = deal([]);
                if isa(state.(name), 'double')
                    if any(strcmpi(name, model.dgVariables))
                        % Get dofs
                        dof = model.getProp(state, [name, 'dof']);
                        % Evaluate at cell cubature points
                        cellValue = model.discretization.evaluateProp(state, dof, 'cell');
                        % Get cell mean
                        cellMean = model.discretization.getCellMean(state, dof);
                        % Evaluate at face cubature points
                        faceValue = model.discretization.evaluateProp(state, dof, 'face');
                    elseif size(double(state.(name)), 1) == model.G.cells.num
                        [fn , index] = model.getVariableField(name, false);
                        if isempty(fn)
                            continue
                        end
                        % Get values
                        v = state.(fn)(:, index);
                        % Repeat to match number cell cubature points
                        cellValue = v(cells,:);
                        % Get cell mean
                        cellMean = v;
                        % Repeat to match number of face cubature points
                        faceValue = v(fcells,:);
                    end
                    if ~isempty(cellValue)
                        % Assign values to the respecive states
                        cellStateDG = model.setProp(cellStateDG, name, cellValue);
                        wellStateDG = model.setProp(wellStateDG, name, cellMean);
                        faceStateDG = model.setProp(faceStateDG, name, faceValue);
                    end
                end
            end
            % Set total saturation
            cellStateDG.sT = getTotalSaturation(cellStateDG.s);
            wellStateDG.sT = getTotalSaturation(wellStateDG.s);
            faceStateDG.sT = getTotalSaturation(faceStateDG.s);
            % Set flag (compositional models)
            if isfield(faceStateDG, 'flag')
                faceStateDG.flag = faceStateDG.flag(fcells);
            end
            % Store cell/well/face states to state
            state.cellStateDG = cellStateDG;
            state.wellStateDG = wellStateDG;
            state.faceStateDG = faceStateDG;
            
        end
        
        function state = assignBaseVariables(model, state)
            
%             names = {'s', 'rs', 'rv'};
            names = fieldnames(state)';
            for name = names
%                 if isfield(state, name{1}) && isfield(state, [name{1}, 'dof'])
                if isfield(state, [name{1}, 'dof'])
                    dof = model.getProp(state, [name{1}, 'dof']);
                    v   = model.discretization.getCellMean(state, dof);
                    state.(name{1}) = v;
                end
            end
             
            if strcmpi(model.formulation, 'totalSaturation')
                if isfield(state, 'sT') && isfield(state, 'sTdof')
                    dof = model.getProp(state, 'stdof');
                    v   = model.discretization.getCellMean(state, dof);
                    state.sT = v;
                end
            end
             
        end
        
        function model = validateModel(model, varargin)
            model = validateModel@TransportModel(model, varargin{:});
                        
            model.parentModel.FluxDiscretization = FluxDiscretizationDG(model.parentModel);
            fp = model.parentModel.FlowPropertyFunctions;
            pvt = fp.getRegionPVT(model.parentModel);
            fp = fp.setStateFunction('PoreVolume', MultipliedPoreVolumeDG(model.parentModel, pvt));
            fp = fp.setStateFunction('GravityPermeabilityGradient', GravityPermeabilityGradientDG(model.parentModel));
            model.parentModel.FlowPropertyFunctions = fp;
            
        end
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            state0 = model.evaluateBaseVariables(state0);
            pmodel = model.parentModel;
            [acc, flux, cellflux, names, types] = pmodel.FluxDiscretization.componentConservationEquations(pmodel, state, state0, dt);
            state.wellStateDG = rmfield(state.wellStateDG, 'FlowProps');
            state.wellStateDG = rmfield(state.wellStateDG, 'FluxProps');
            src = pmodel.FacilityModel.getComponentSources(state.wellStateDG);
            % Treat source or bc terms
            if ~isempty(drivingForces.bc) || ~isempty(drivingForces.src)
                fluxBC  = model.computeBoundaryConditions(state, state0, dt, drivingForces.bc);
            end
            % Assemble equations and add in sources
            if strcmpi(model.formulation, 'missingPhase')
                % Skip the last phase! Only mass-conservative for
                % incompressible problems
                acc   = acc(1:end-1);
                flux  = flux(1:end-1);
                names = names(1:end-1);
                types = types(1:end-1);
            end
            d        = model.discretization;
            d.nDof   = state.nDof;
            d.dofPos = state.dofPos;
            psi      = d.basis.psi;
            gradPsi  = d.basis.gradPsi;
            ixw      = d.getDofIx(state, 1, src.cells);
            ix       = d.getDofIx(state, Inf);
            d.sample = state.cellStateDG.s{1}(ix)*0;
            eqs      = cell(1, numel(acc));
            state.wellStateDG.cells = (1:pmodel.G.cells.num)';
            
            cells  = rldecode((1:pmodel.G.cells.num)', state.nDof, 1);
            pv     = pmodel.operators.pv(cells);
            rhoS   = pmodel.getSurfaceDensities();
            cnames = pmodel.getComponentNames();
            
            for i = 1:numel(acc)
                eqs{i} = d.inner(acc{i}     , psi    , 'dV') ...
                       - d.inner(cellflux{i}, gradPsi, 'dV') ...
                       + d.inner(flux{i}    , psi    , 'dS');
                if ~isempty(drivingForces.bc)
                    eqs{i} = eqs{i} + d.inner(fluxBC{i}, psi, 'dS', drivingForces.bc.face);
                end
                if ~isempty(src.cells)
                    eqs{i}(ixw) = eqs{i}(ixw) - src.value{i};
                end
                if ~pmodel.useCNVConvergence
                    sub = strcmpi(names{i}, cnames);
                    eqs{i} = eqs{i}.*(dt./(pv.*rhoS(sub)));
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function q = computeBoundaryConditions(model, state, state0, dt, bc)
            
            bcState = model.parentModel.FluxDiscretization.buildFlowState(model, state, state0, dt);
            faces = bc.face;
            [~, x, ~, fNo] = model.discretization.getCubature(faces, 'face');
            cNo = sum(model.parentModel.G.faces.neighbors(fNo,:),2);
            names = fieldnames(bcState);
            for k = 1:numel(names)
                name = names{k};
                if numel(name) > 3 && strcmp(name(end-2:end), 'dof')
                    % Get dofs
                    dof = model.getProp(bcState, name);
                    if iscell(dof)
                        v = cell(1,numel(dof));
                        for i = 1:numel(dof)
                            v{i} = model.discretization.evaluateDGVariable(x, cNo, state, dof{i});
                        end
                    else
                        v = model.discretization.evaluateDGVariable(x, cNo, state, dof);
                    end
                    n = name(1:end-3);
                    % Evaluate at boundary face cubature points
                    bcState = model.setProp(bcState, n, v);
                end
            end
            bcState.cells = sum(model.G.faces.neighbors(fNo,:), 2);
            bcState.faces = fNo;
            
            q = computeBoundaryFluxesDG(model.parentModel, bcState, bc);
            
        end
        
        %-----------------------------------------------------------------%
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            [model, state] = prepareTimestep@TransportModel(model, state, state0, dt, drivingForces);
            state = assignDofFromState(model.discretization, state, {'pressure'});
        end
        
        %-----------------------------------------------------------------%
        function [restVars, satVars, wellVars] = splitPrimaryVariables(model, vars)
            vars = cellfun(@(n) n(1:end-3), vars, 'UniformOutput', false);
            [restVars, satVars, wellVars] = model.parentModel.splitPrimaryVariables(vars);
            restVars = cellfun(@(n) [n, 'dof'], restVars, 'UniformOutput', false);
            satVars = cellfun(@(n) [n, 'dof'], satVars, 'UniformOutput', false);
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
                        
            % Remove DG states
            state = rmfield(state, 'cellStateDG');
            state = rmfield(state, 'faceStateDG');
            state = rmfield(state, 'wellStateDG');
            
            if 0% strcmpi(class(model.parentModel), 'GenericBlackOilModel') ...
%                  && model.parentModel.disgas || model.parentModel.vapoil
%                 [state, report] = model.updateStateBO(state, problem, dx, drivingForces);
            else
                state0 = state;
                [restVars, satVars] = model.splitPrimaryVariables(problem.primaryVariables);
                % Update saturation dofs
                state = model.updateSaturations(state, dx, problem, satVars);
                % Update non-saturation dofs
                state = model.updateDofs(state, dx, problem, restVars);
                % Update cell averages from dofs
                state  = model.assignBaseVariables(state);
                report = [];

                if 1
                % Compute dx for cell averages
                dx0 = model.getMeanIncrement(state, state0, problem.primaryVariables);
                % Let parent model do its thing
                problem0 = problem;
                problem0.primaryVariables = cellfun(@(n) n(1:end-3), problem0.primaryVariables, 'UniformOutput', false);
                [state0_corr, report] = updateState@TransportModel(model, state0, problem0, dx0, drivingForces);
                % Correct updates in dofs according to parent model
                dx0_corr = model.getMeanIncrement(state0_corr, state0, problem.primaryVariables);
                cells    = rldecode((1:model.G.cells.num)', state.nDof, 1);
                frac     = cellfun(@(x,y) x(cells)./y(cells), dx0_corr, dx0, 'UniformOutput', false);
                for i = 1:numel(frac)
                    frac{i}(~isfinite((frac{i}))) = 1;
                end
                dx_corr  = cellfun(@(dx, f) dx.*f, dx, frac, 'UniformOutput', false);
                % Update saturation dofs
                state = model.updateSaturations(state0, dx_corr, problem, satVars);
                % Update non-saturation dofs
                state = model.updateDofs(state, dx_corr, problem, restVars);
                % Compositional models
                if isfield(state0_corr, 'eos')
                    state.eos = state0_corr.eos;
                    varsEOS = {'Kdof', 'Ldof', 'xdof', 'ydof', 'Z_Ldof', 'Z_Vdof', 'sodof', 'sgdof'};
                    dxEOS   = model.getMeanIncrement(state0_corr, state, varsEOS);
                    for i = 1:numel(varsEOS)
                        [fn, index] = model.getVariableField(varsEOS{i});
                        ix = model.discretization.getDofIx(state, 1);
                        state.(fn)(ix,index) = state.(fn)(ix,index) + dxEOS{i};
                    end
%                     state   = model.updateDofs(state, dxEOS, problem, varsEOS);
                end
                % Update cell averages from dofs
                state = model.assignBaseVariables(state);
                end
            
            end
            
        end
        
        % --------------------------------------------------------------------%
        function [state, report] = updateStateBO(model, state, problem, dx, drivingForces)
            vars    = problem.primaryVariables;
            removed = false(size(vars));
            pmodel  = model.parentModel;
            cells   = rldecode((1:model.G.cells.num)', state.nDof, 1);
            if pmodel.disgas || pmodel.vapoil
                % The VO model is a bit complicated, handle this part
                % explicitly.
                state0 = state;
                state = model.initStateFunctionContainers(state);

                state = pmodel.updateStateFromIncrement(state, dx, problem, 'pressure', pmodel.dpMaxRel, pmodel.dpMaxAbs);
                state = pmodel.capProperty(state, 'pressure', pmodel.minimumPressure, pmodel.maximumPressure);

                [vars, ix] = model.stripVars(vars, 'pressure');
                removed(~removed) = removed(~removed) | ix;

                % Black oil with dissolution
                [so, sg] = model.getProps(state, 'so', 'sg');
                if pmodel.water
                    sw  = model.getProp(state, 'sw');
                    dsw = model.getIncrement(dx, problem, 'swdof');
                else
                    sw  = 0;
                    dsw = 0;
                end
                % Magic status flag, see inside for doc
                st0 = pmodel.getCellStatusVO(state0, so, sw, sg);
                st = st0;
                for j = 1:numel(st)
                    if numel(st{j}) == model.G.cells.num
                        st{j} = st{j}(cells);
                    end
                end
                dr = model.getIncrement(dx, problem, 'xdof');
                % Interpretation of "gas" phase varies from cell to cell, remove
                % everything that isn't sG updates
                dsg = st{3}.*dr - st{2}.*dsw;

                drsMaxAbs = pmodel.drsMaxAbs/model.discretization.basis.nDof;
                if pmodel.disgas
                    rsMax = pmodel.getProp(state, 'rsMax');
                    rsMax = rsMax(cells);
                    drs_rel = rsMax.*pmodel.drsMaxRel/model.discretization.basis.nDof;
                    drs = min(drsMaxAbs, drs_rel);
                    state = model.updateStateFromIncrement(state, st{1}.*dr, problem, ...
                                                           'rsdof', inf, drs);
                    state.rs = model.discretization.getCellMean(state, state.rsdof);
                end
                
                if 0
                    rs    = model.getProp(state, 'rs');
                    rsSat = model.parentModel.getProp(state, 'rsMax');
                    rsdof = model.getProp(state, 'rsdof');
                    [rsMin, rsMax] = model.discretization.getMinMax(state, rsdof);
                    rsMin(rs > rsSat) = rsSat(rs > rsSat);
                    rsMax(rs < rsSat) = rsSat(rs < rsSat);
                    state = model.discretization.limiter{2}(state, 'rs', [rsMin, rsMax]);
                end

                if pmodel.vapoil
                    rvMax = pmodel.getProp(state, 'rvMax');
                    drv_rel = rvMax.*pmodel.drsMaxRel/model.discretization.basis.nDof;
                    drs = min(drsMaxAbs, drv_rel);
                    state = model.updateStateFromIncrement(state, st{2}.*dr, problem, ...
                                                           'rvdof', inf, drs);
                end

                dso = -(dsg + dsw);
                nPh = nnz(pmodel.getActivePhases());

                ds = zeros(numel(dso), nPh);
                phIndices = pmodel.getPhaseIndices();
                if pmodel.water
                    ds(:, phIndices(1)) = dsw;
                end
                if pmodel.oil
                    ds(:, phIndices(2)) = dso;
                end
                if pmodel.gas
                    ds(:, phIndices(3)) = dsg;
                end
                
                dsMaxAbs = pmodel.dsMaxAbs/model.discretization.basis.nDof;
                state = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, dsMaxAbs);
                state = model.assignBaseVariables(state);
                
                kr = pmodel.FlowPropertyFunctions.RelativePermeability;
                state = kr.applyImmobileChop(model, state, state0);

                % We should *NOT* be solving for oil saturation for this to make sense
                assert(~any(strcmpi(vars, 'sodof')));
                
%                 problem0 = problem;
                names = {'swdof', 'sodof', 'sgdof'};
                if pmodel.disgas
                    names{end+1} = 'rsdof';
                end
                if pmodel.vapoil
                    names{end+1} = 'rvdof';
                end
                dx0 = model.getMeanIncrement(state, state0, names);
                
                state_corr = computeFlashBlackOil(state, state0, pmodel, st0);
                state_corr.s = bsxfun(@rdivide, state_corr.s, sum(state_corr.s, 2));

                problem0 = problem;
                problem0.primaryVariables = names;
                dx0_corr = model.getMeanIncrement(state_corr, state0, names);
                
                frac = cell(1, numel(dx0));
                [dx_corr, dx] = deal(cell(numel(names),1));
                for i = 1:numel(names)
                    v0 = model.getProp(state0, names{i});
                    v  = model.getProp(state , names{i});
                    dx{i} = v - v0;
                end
                for i = 1:numel(frac)
                    if  ~isempty(dx0{i})
                        f = dx0_corr{i}(cells)./dx0{i}(cells);
                        f(~isfinite(f)) = 1;
                        dx_corr{i} = dx{i}.*f;
                    end
                end
                % Update saturation dofs
                state = model.updateDofs(state0, dx_corr, problem0, names);
                state.status = state_corr.status;
                % Update cell averages from dofs
                state = model.assignBaseVariables(state);
                
                if 0
                    rsSat = model.parentModel.getProp(state, 'RsMax');
                    rs    = model.getProp(state, 'rs');
                    rsdof = model.getProp(state, 'rsdof');
                    [rsMin, rsMax] = model.discretization.getMinMax(state, rsdof);
                    rsMin(rs > rsSat) = rsSat(rs > rsSat);
                    rsMax(rs < rsSat) = rsSat(rs < rsSat);
                    state = model.discretization.limiter{2}(state, 'rs', [rsMin, rsMax]);
                end
                
                %  We have explicitly dealt with rs/rv properties, remove from list
                %  meant for autoupdate.
                [vars, ix] = model.stripVars(vars, {'swdof', 'sodof', 'sgdof', 'rsdof', 'rvdof', 'xdof'});
                removed(~removed) = removed(~removed) | ix;
            end

            % We may have solved for a bunch of variables already if we had
            % disgas / vapoil enabled, so we remove these from the
            % increment and the linearized problem before passing them onto
            % the generic reservoir update function.
            problem.primaryVariables = vars;
            dx(removed) = [];
            report = [];
%             [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
        end
        
        %-----------------------------------------------------------------%
        function dx = getMeanIncrement(model, state, state0, vars)
            
%             vars = problem.primaryVariables;
            dx   = cell(numel(vars),1);
            for i = 1:numel(vars)
                vn = vars{i}(1:end-3);
                v  = model.getProp(state, vn);
                v0 = model.getProp(state0, vn);
                dx{i} = v - v0;
            end

        end
        
        % ----------------------------------------------------------------%
        function state = updateDofs(model, state, dx, problem, dofVars, dvMaxAbs)
            
            for i = 1:numel(dofVars)
                dvMaxAbs = inf;
                if strcmpi(dofVars{i}, 'sTdof') && 0
                    dvMaxAbs = 0.2;
                end
                state = updateStateFromIncrement(model, state, dx, problem, dofVars{i}, inf, dvMaxAbs);
            end
            
        end
        
        % ----------------------------------------------------------------%
        function state = updateSaturations(model, state, dx, problem, satVars)

            if nargin < 5
                % Get the saturation names directly from the problem
                [~, satVars] = ...
                    splitPrimaryVariables(model, problem.primaryVariables);
            end
            if isempty(satVars)
                % No saturations passed, nothing to do here.
                return
            end
            % Solution variables should be saturations directly, find the
            % missing link
            saturations0 = lower(model.parentModel.getSaturationVarNames);
            saturations  = cellfun(@(n) [n, 'dof'], saturations0, 'uniformOutput', false);
            fillsat = setdiff(saturations, lower(satVars));
            nFill = numel(fillsat);
            assert(nFill == 0 || nFill == 1)
            if nFill == 1
                % Fill component is whichever saturation is assumed to fill
                % up the rest of the pores. This is done by setting that
                % increment equal to the negation of all others so that
                % sum(s) == 0 at end of update
                fillsat = fillsat{1};
                solvedFor = ~strcmpi(saturations, fillsat);
            else
                % All saturations are primary variables. Sum of saturations is
                % assumed to be enforced from the equation setup
                solvedFor = true(numel(saturations), 1);
            end
            ds = zeros(sum(state.nDof), numel(saturations));
            
            tmp = 0;
            ix = model.discretization.getDofIx(state, Inf);
            for phNo = 1:numel(saturations)
                if solvedFor(phNo)
                    v = model.getIncrement(dx, problem, saturations{phNo});
                    ds(ix, phNo) = v;
                    if nFill > 0
                        % Saturations added for active variables must be subtracted
                        % from the last phase
                        tmp = tmp - v;
                    end
                end
            end
            ds(ix, ~solvedFor) = tmp;
            % We update all saturations simultanously, since this does not bias the
            % increment towards one phase in particular.
            if 1
                dsAbsMax = model.parentModel.dsMaxAbs/model.discretization.basis.nDof;
            else
                dsAbsMax = model.parentModel.dsMaxAbs/min(model.discretization.basis.nDof, 3);
            end
            state = model.updateStateFromIncrement(state, ds, problem, 'sdof', Inf, dsAbsMax);
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            state.FacilityFluxProps = state.wellStateDG.FacilityFluxProps;
            [state, report] = updateAfterConvergence@TransportModel(model, state0, state, dt, drivingForces);
            state = rmfield(state, 'cellStateDG');
            state = rmfield(state, 'faceStateDG');
            state = rmfield(state, 'wellStateDG');
            
            propfn = model.parentModel.getStateFunctionGroupings();
            d = model.discretization;
            d.nDof = state.nDof;
            d.dofPos = state.dofPos;
            ix = d.getDofIx(state, 1, Inf);
            psi    = model.discretization.basis.psi(1);
            d.sample = state.sdof(:,1);
            for i = 1:numel(propfn)
                p = propfn{i};
                struct_name = p.getStateFunctionContainerName();
                names = p.getNamesOfStateFunctions();
                if isfield(state, struct_name)
                    for j = 1:numel(names)
                        name = names{j};
                        if ~isempty(state.(struct_name).(name))
                            v = state.(struct_name).(name);
                            nph = numel(v);
                            for ph = 1:nph
                                v{ph} = d.inner(v{ph}, psi, 'dV');
                                v{ph} = v{ph}(ix);
                            end
                            state.(struct_name).(name) = v;
                        end
                    end
                end
            end
            
            if ~isempty(model.limiters)
                
                if model.storeUnlimited
                    state.ul = state;
                    if isfield(state.ul, 'ul')
                        state.ul = rmfield(state.ul, 'ul');
                    end
                end
                
                for l = 1:numel(model.limiters)
                    limiter = model.limiters(l);
                    for v = 1:numel(limiter.variables)
                        state = limiter.function(state, limiter.variables{v}, limiter.tol, limiter.limits{v});
                    end
                end
                
                if 0
                    if isa(model.parentModel, 'GenericBlackOilModel')
                        if model.parentModel.disgas
                            rsSat = model.parentModel.getProp(state, 'RsMax');
                            rs    = model.getProp(state, 'rs');
                            rsdof = model.getProp(state, 'rsdof');
                            [rsMin, rsMax] = model.discretization.getMinMax(state, rsdof);
                            rsMin(rs > rsSat) = rsSat(rs > rsSat);
                            rsMax(rs < rsSat) = rsSat(rs < rsSat);
                            state = model.discretization.limiter{2}(state, 'rs', [rsMin, rsMax]);
                            state = model.discretization.limiter{1}(state, 'rs');
                        end
                        if model.parentModel.vapoil
                            state = model.discretization.limiter{1}(state, 'rv');
                        end
                    end
                end
            end

        end
        
    end
    
end

function sT = getTotalSaturation(s)
    if iscell(s)
        sT  = 0;
        nph = numel(s);
        for i = 1:nph
            sT = sT + s{i};
        end
    else
        sT = sum(s,2);
    end
end