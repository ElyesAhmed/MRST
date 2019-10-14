classdef SpatialDiscretization
        
    properties
        G
        internalConn
        N
    end
    
    methods 
        function disc = SpatialDiscretization(G, varargin)
            disc.G            = G;
            N                 = G.faces.neighbors;
            disc.internalConn = all(N ~= 0, 2);
            disc.N            = N(disc.internalConn,:);
        end
        
    end
    
end