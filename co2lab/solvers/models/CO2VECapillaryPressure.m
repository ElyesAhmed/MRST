classdef CO2VECapillaryPressure < StateFunction
   
    properties
    end
    
    methods
        function prop = CO2VECapillaryPressure(model, varargin)
            prop = prop@StateFunction(model, varargin{:});
            prop = prop.dependsOn({'s', 'sGmax', 'pressure'}, 'state');
            prop.label = 'p_{c}';
        end
        
        function pc = evaluateOnDomain(prop, model, state)
            sG = model.getProp(state, 'sg');
            sGmax = model.getProp(state, 'sGmax');
            pressure = model.getProp(state, 'pressure');
            pc = prop.evaluateFluid(model, 'pcWG', sG, pressure, 'sGmax', sGmax);
            pc = {0*pc, pc}; % we need one column for water, another for CO2,
                             % indicating what to be added to the reference 
                             % pressure to get the phase pressure.  Since
                             % reference pressure is the same as water
                             % pressure here, the corresponding column is zero.
        end
        
        function res = pcPresent(prop, model)
            res = true;
        end
    end
end
