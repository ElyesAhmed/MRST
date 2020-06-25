classdef MechModel < PhysicalModel

    properties
        
        mech;
        % Structure that contains the mechanical parameter of the system
        % This structure should contain the fields prop
        % prop :  - lambda, first Lame parameter
        %         - mu    , second Lame parameter
        %
        % and the loading structure loadstruct
        % loadstruct:  - bc      , Dirichlet boundary conditions 
        %              - extforce, Neumann boundary conditions (external forces)
        %              - force   , Volumetric boundary conditions

        MechPropertyFunctions;
        
    end

    methods
        
        function model = MechModel(G, mech, varargin)
            
            model = model@PhysicalModel(G, varargin{:});

            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end

            % Physical properties of rock and fluid
            model.mech  = mech;

            model.MechPropertyFunctions = MechPropertyFunctions(model);
            
            model.operators = setupMpsaOperators(model);
            
        end

        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            [eqs, names, types, state] = mechEquations(model, state0, state, dt, drivingForces, varargin);
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            [u, lambdamech] = model.getProps(state, 'displacement', 'lambdamech');
            vars = {u, lambdamech};
            names = {'displacement', 'lambdamech'};
            origin = {class(model)};
        end
        
        
        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@PhysicalModel(model);
            containers = [containers, {model.MechPropertyFunctions}];
        end
        
        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@PhysicalModel(model, varargin{:});
            model.MechPropertyFunctions = MechPropertyFunctions(model);
        end
    
        function u = vectorDisplacement(model, state)
            u = model.getProp(state, 'displacement');
            dim = model.G.griddim;
            u = reshape(u, dim, [])';
        end

        function state = solveMechanics(model)
            
            state.u = model.operators.rhs{1}; % Dummy values, just used to get the correct dimension.
            state.lambdamech = model.operators.rhs{2}; % Dummy values, just used to get the correct dimension
            
            solver = NonLinearSolver(); 
            drivingForces = [];
            [state, failure, report] = solveMinistep(solver, model, state, state, 0, drivingForces);
            
        end
    
        function [fn, index] = getVariableField(model, name, varargin)
        % Get the index/name mapping for the model (such as where pressure or water saturation is located in state)
            switch(lower(name))
              case {'u', 'displacement'}
                % Displacement
                fn = 'u';
                index = ':';
              case {'lambdamech'}
                % Lagrangian variables for corresponding to the boundary conditions
                fn = 'lambdamech';
                index = ':';
              otherwise
                % This will throw an error for us
                [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
            end
        end

    end
end
