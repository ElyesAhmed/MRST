classdef FixedTotalVelocityDG < FixedTotalFluxDG
    
    properties
    end
    
    methods
        function gp = FixedTotalVelocityDG(model, varargin)
            gp@FixedTotalFluxDG(model, varargin{:});
        end
        function vT = evaluateOnDomain(prop, model, state)
            vT = model.disc.velocityInterp.faceFlux2cellVelocity(sum(state.(prop.fluxfield),2));
            vT = vT(state.cells,:);
        end
    end
end