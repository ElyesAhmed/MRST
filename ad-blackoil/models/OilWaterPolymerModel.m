classdef OilWaterPolymerModel < TwoPhaseOilWaterModel
    % Oil/water/polymer system
    % This model is a two phase oil/water model, extended with the polymer
    % phase in addition.
    
    properties
        % Polymer present
        polymer
        
    end
    
    methods
        function model = OilWaterPolymerModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            
            % This is the model parameters for oil/water/polymer
            model.polymer = true;
            
            model.wellVarNames = {'qWs', 'qOs', 'qWPoly', 'bhp'};
            
            model = merge_options(model, varargin{:});
            
        end
        
        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterPolymer(state0, state, ...
                model, dt, drivingForces, varargin{:});
        end
        
        function [state, report] = updateState(model, state, problem, ...
                dx, drivingForces)
            [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
               state, problem,  dx, drivingForces);
            
            if model.polymer
                c = model.getProp(state, 'polymer');
                c = min(c, model.fluid.cmax);
                state = model.setProp(state, 'polymer', max(c, 0) );
            end
        end
        
        function state = updateAfterConvergence(model, state0, state, dt, drivingForces)
            state = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, dt, drivingForces);
            if model.polymer
                c     = model.getProp(state, 'polymer');
                cmax  = model.getProp(state, 'polymermax');
                state = model.setProp(state, 'polymermax', max(cmax, c));
            end
        end

        
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case {'polymer'}
                    index = 1;
                    fn = 'c';
                case {'polymermax'}
                    index = 1;
                    fn = 'cmax';
                otherwise
                    [fn, index] = getVariableField@ReservoirModel(...
                                    model, name);
            end
        end
        
        
    end
end
