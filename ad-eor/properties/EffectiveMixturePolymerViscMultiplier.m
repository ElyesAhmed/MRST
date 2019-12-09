classdef EffectiveMixturePolymerViscMultiplier < StateFunction
    properties
    end

    methods
        function gp = EffectiveMixturePolymerViscMultiplier(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'polymer'}, 'state'); % check mechanism
        end

        function muWeffMult = evaluateOnDomain(prop, model, state)
            cp   = model.getProp(state, 'polymer');            
            fluid = model.fluid;
            mixpar = fluid.mixPar;
            cpbar   = cp/fluid.cpmax;
            a = fluid.muWMult(fluid.cpmax).^(1-mixpar);
            b = 1./(1 - cpbar + cpbar./a);
            % The viscosity multiplier only results from the polymer mixing.
            muWeffMult = b.*fluid.muWMult(cp).^mixpar;            
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
