classdef PhaseMixingCoefficientsLV < StateFunction
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = PhaseMixingCoefficientsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'temperature', 'x', 'y'}, 'state');
        end

        function v = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            [p, T, x, y] = model.getProps(state,...
                'pressure', 'temperature', 'liquidMoleFractions', 'vaporMoleFractions');
            acf = eos.fluid.acentricFactors;
            
            [A_ij, Bi] = eos.getMixingParameters(p, T, acf, iscell(x));
            [Si_L, A_L, B_L] = eos.getPhaseMixCoefficients(x, A_ij, Bi);
            twoPhase = model.getTwoPhaseFlag(state);
            if prop.useCompactEvaluation && iscell(y) && ~all(twoPhase)
                Si_V = Si_L;
                A_V = A_L;
                B_V = B_L;
                if any(twoPhase)
                    A_ij_2ph = cellfun(@(x) x(twoPhase), A_ij, 'UniformOutput', false);
                    Bi_2ph = cellfun(@(x) x(twoPhase), Bi, 'UniformOutput', false);
                    y_2ph = cellfun(@(x) x(twoPhase), y, 'UniformOutput', false);
                    [Si_V_2ph, A_V(twoPhase), B_V(twoPhase)] = eos.getPhaseMixCoefficients(y_2ph, A_ij_2ph, Bi_2ph);
                    for i = 1:numel(Si_V_2ph)
                        Si_V{i}(twoPhase) = Si_V_2ph{i};
                    end
                end
            else
                [Si_V, A_V, B_V] = eos.getPhaseMixCoefficients(y, A_ij, Bi);
            end
            
            L = struct('Si', {Si_L}, 'A', {A_L}, 'B', {B_L}, 'Bi', {Bi}, 'A_ij', {A_ij});
            V = struct('Si', {Si_V}, 'A', {A_V}, 'B', {B_V}, 'Bi', {Bi}, 'A_ij', {A_ij});
            if model.water
                v = {struct(); L; V};
            else
                v = {L; V};
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
