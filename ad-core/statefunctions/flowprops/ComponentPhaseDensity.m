classdef ComponentPhaseDensity < StateFunction & ComponentProperty
    % Component density in each cell for each phase
    properties
        includeStandard = true;
    end
    
    methods
        function gp = ComponentPhaseDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp@ComponentProperty(model, 'getComponentDensity');
            gp.label = '\rho_\alpha x_{i,\alpha}';
            if gp.includeStandard
                gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
            end
        end
        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            v = cell(ncomp, nph);
            extra = prop.getExtraArguments(model, state);
            for c = 1:ncomp
                v(c, :) = model.Components{c}.getComponentDensity(model, state, extra{:});
            end
        end
        
        function extra = getExtraArguments(prop, model, state)
            if prop.includeStandard
                [rho] = prop.getEvaluatedExternals(model, state, 'Density');
                extra = {rho};
            else
                extra = {};
            end
        end
    end
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
