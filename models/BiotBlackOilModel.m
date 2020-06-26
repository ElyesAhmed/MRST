classdef BiotBlackOilModel < GenericBlackOilModel

    properties
        
        mech

        eta
        bcetazero
        
        BiotPropertyFunctions
        MechPropertyFunctions        
        
    end

    methods
        
        function model = BiotBlackOilModel(G, rock, fluid, mech, varargin)
            
            model = model@GenericBlackOilModel(G, rock, fluid, varargin{:});
            model = merge_options(model, varargin{:});
            
            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end
            % Physical properties of rock and fluid
            model.mech  = mech;
        
            model.eta = 0;
            model.bcetazero = false;
            
            % Add mechanical operators  
            model.operators = setupBiotAdOperators(model);

            model.BiotPropertyFunctions = BiotPropertyFunctions(model);
            model.MechPropertyFunctions = MechPropertyFunctions(model);
            
        end
        
        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@GenericBlackOilModel(model);
            containers = [containers, {model.BiotPropertyFunctions, model.MechPropertyFunctions}];
        end

        function model = setupStateFunctionGroupings(model, varargin) 
            
            model.PVTPropertyFunctions = []; % make sure this ir reset
            model = setupStateFunctionGroupings@GenericBlackOilModel(model, varargin{:});
           
            fluxprops = model.FluxDiscretization;
            biotprops = model.BiotPropertyFunctions; 
            pvtprops  = model.PVTPropertyFunctions; 
            mprops  = model.MechPropertyFunctions;
            
            pv = pvtprops.getStateFunction('PoreVolume');
            
            biotprops = biotprops.setStateFunction('BasePoreVolume'      , pv);
            biotprops = biotprops.setStateFunction('Dilatation'          , BiotBlackOilDilatation(model));
            pvtprops  =  pvtprops.setStateFunction('PoreVolume'          , BiotPoreVolume(model));
            mprops    =    mprops.setStateFunction('FaceNodeDisplacement', BiotFaceNodeDisplacement(model));
            fluxprops = fluxprops.setStateFunction('PermeabilityPotentialGradient', MpfaKgrad(model));
            fluxprops = fluxprops.setStateFunction('PhaseUpwindFlag', MpfaPhaseUpwindFlag(model));
            
            model.BiotPropertyFunctions = biotprops;
            model.PVTPropertyFunctions  = pvtprops;
            model.MechPropertyFunctions = mprops;
            model.FluxDiscretization    = fluxprops;
            
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            [eqs, names, types, state] = biotBlackOilEquations(model, state0, state, dt, drivingForces, varargin);
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)

            [vars, names, origin] = getPrimaryVariables@GenericBlackOilModel(model, state);
            [u, bp, lm] = model.getProps(state, 'u', 'biotpressure', 'lambdamech');
            vars = [vars, {u, bp, lm}];
            names = [names, {'displacement', 'biotpressure', 'lambdamech'}];

            origin = [origin, {class(model), class(model), class(model)}];
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
            
            switch(lower(name))
              case {'u', 'displacement'}
                % Displacement
                fn = 'u';
                index = ':';
              case {'lambdamech'}
                % Lagrangian variables for corresponding to the boundary conditions
                fn = lower(name);
                index = ':';
              case {'biotpressure', 'bp'}
                fn = 'biotpressure';
                index = ':';
              otherwise
                [fn, index] = getVariableField@GenericBlackOilModel(model, name, varargin{:});
            end

        end

    end

end
