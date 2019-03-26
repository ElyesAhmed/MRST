classdef MultipliedPoreVolume < GridProperty
    properties
    end
    
    methods
        function gp = MultipliedPoreVolume(model, varargin)
            gp@GridProperty(model, varargin{:});
            if isfield(model.fluid, 'pvMultR')
                gp = gp.dependsOn({'pressure'}, 'state');
            end
        end
        function pv = evaluateOnDomain(prop, model, state)
            f = model.fluid;
            pv = model.operators.pv;
            if isfield(f, 'pvMultR')
                p = model.getProp(state, 'pressure');
                pvMult = prop.evaluateFunctionOnGrid(f.pvMultR, p);
                pv = pv.*pvMult;
            end
        end
    end
end