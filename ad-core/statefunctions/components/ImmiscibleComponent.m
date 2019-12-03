classdef ImmiscibleComponent < ComponentImplementation
    % Specialized interface for immiscible component
    properties
        phaseIndex % Index of phase this component belongs to
    end
    
    methods
        function c = ImmiscibleComponent(name, phase)
            c@ComponentImplementation(name);
            c.phaseIndex = phase;
            c = c.dependsOn('Density');
        end
        
        function c = getComponentDensity(component, model, state, varargin)
            c = getComponentDensity@ComponentImplementation(component, model, state, varargin{:});
            rho = model.getProp(state, 'Density');
            c{component.phaseIndex} = rho{component.phaseIndex};
        end
        
        function c = getPhaseComposition(component, model, state, varargin)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            c{component.phaseIndex} = 1;
        end
        
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            c = component.getPhaseComposition(model, state);
        end
        
        function c = getPhaseComponentFractionWell(component, model, state, W)
            % Get the fraction of the component in each phase (when
            % injecting from outside the domain)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            comp_i = vertcat(W.compi);
            index = component.phaseIndex;
            ci = comp_i(:, index);
            if any(ci ~= 0)
                c{index} = ci;
            end
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
