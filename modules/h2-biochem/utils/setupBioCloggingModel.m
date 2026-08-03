function [model, poro0, perm0] = setupBioCloggingModel(model, nbact0, nc, cp, clogModel)
% setupBioCloggingModelMulti -- Add multi-species bio-clogging effects
%
% PARAMETERS:
%   nbact0 - Cell array of initial bacterial concentrations (e.g., {nbact1_0, nbact2_0})
%   nc     - Vector of characteristic concentrations matching the cell array length [nc1, nc2]
%   cp     - Vector of clogging strengths matching the cell array length [cp1, cp2]

if nargin < 5
    clogModel = true;
end

poro0 = model.rock.poro;
perm0 = model.rock.perm(:, 1);

if clogModel
    % 1. Compute the cumulative initial "scale" factor across all species
    if ~iscell(nbact0)
        if numel(nc) == 1
            nbact0 = {nbact0};
        elseif isvector(nbact0) && numel(nbact0) == numel(nc)
            nbact0 = num2cell(nbact0);
        else
            assert(size(nbact0, 2) == numel(nc), ...
                'nbact0 must have one column per bacterial species.');
            nbact0 = num2cell(nbact0, 1);
        end
    end
    num_species = numel(nbact0);
    assert(numel(nc) == num_species && numel(cp) == num_species, ...
        'nc and cp must contain one value per bacterial species.');
    scale_sum = 0;
    for i = 1:num_species
        scale_sum = scale_sum + cp(i) * (nbact0{i} ./ nc(i)).^2;
    end
    scale = 1 + scale_sum;

    % 2. Define the dynamic pore volume multiplier for multiple species
    pvMult_nbact = @(varargin) 1 ./ (1 + scale .* evalCumulativeClog(varargin, nc));

    % 3. Assign porosity handles
    model.fluid.pvMultR = @(p, varargin) pvMult_nbact(varargin{:});

    poroFun = @(p, varargin) poro0 .* pvMult_nbact(varargin{:});
    model.rock.poro = poroFun;

    % 4. Define permeability update function (Kozeny–Carman)
    tauFun = @(p, varargin) ((1 - poro0) ./ (1 - poroFun(p, varargin{:}))).^2 .* ...
        (poroFun(p, varargin{:}) ./ poro0).^3;
    permFun = @(p, varargin) perm0 .* tauFun(p, varargin{:});
    model.rock.perm = permFun;
else
    model.rock.poro = poro0;
    model.rock.perm = perm0;
    model.fluid.pvMultR = @(varargin) 1;
end

% Rock handles determine which dynamic state functions are registered.
% Rebuild operators and state functions after replacing those handles.
model = model.setupOperators();
model.FlowDiscretization = BiochemicalFlowDiscretization(model);
model = model.setupStateFunctionGroupings();

end

% Helper function to loop over cell contents during simulation solver steps
function total_clog = evalCumulativeClog(nbact_cell, nc_vec)
total_clog = 0;
for i = 1:numel(nbact_cell)
    total_clog = total_clog + (nbact_cell{i} ./ nc_vec(i)).^2;
end
end
