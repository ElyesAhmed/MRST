function v = discreteDivergence(acc, N, v, nc, nf, sortIx, C, prelim, useMex)
% Discrete divergence for the GenericAD library

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

    hasAcc = not(isempty(acc));
    if isa(v, 'GenericAD')
        v.val = accumulate(N, value(v), nc);
        if hasAcc && isa(acc, 'GenericAD')
            % Both present, both are AD
            for i = 1:numel(v.jac)
                v.jac{i} = accDivJac(acc.jac{i}, v.jac{i}, N, nc, nf, sortIx, C, prelim, useMex);
            end
            v.val = v.val + acc.val;
        else
            for i = 1:numel(v.jac)
                v.jac{i} = divJac(v.jac{i}, N, nc, nf, sortIx, C, prelim, useMex);
            end
            if hasAcc
                % Acc is not AD
                v.val = v.val + acc;
            end
        end
    else
        assert(isnumeric(v), 'Expected numeric vector, but got ''%s''\n', class(v))
        v = accumulate(N, v, nc);
        if hasAcc
            v = v + acc;
        end
    end
end

function v = accumulate(N, v, nc)
    v = accumarray(N(:, 1), v, [nc, 1]) - accumarray(N(:, 2), v, [nc, 1]);
end

function jac = divJac(jac, N, nc, nf, sortIx, C, prelim, useMex)
    if issparse(jac)
        if nnz(jac) > 0
            if isempty(C)
                C  = sparse(N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);
            end
            jac = C*jac;
        else
            jac = sparse([], [], [], nc, matrixDims(jac, 2));
        end
    elseif jac.isZero
            jac = sparse([], [], [], nc, prod(jac.dim));
        return
    else
        if useMex && (isempty(jac.parentSubset) || all(jac.parentSubset == (1:jac.dim(1))'))
            jac = mexDiscreteDivergenceJac([], jac.diagonal, N, prelim.facePos, prelim.faces, prelim.cells, prelim.cellIndex);
        else
            jac = sortIx.C*jac.sparse();
        end
    end
end

function jac = accDivJac(acc, jac, N, nc, nf, sortIx, C, prelim, useMex)
    if issparse(jac)
        if nnz(jac) > 0
            if isempty(C)
                C  = sparse(N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);
            end
            jac = C*jac + acc;
        else
            jac = acc;
        end
    elseif jac.isZero
            jac = acc;
        return
    else
        if useMex && (isempty(jac.parentSubset) || (numel(jac.parentSubset) == jac.dim(1)) && all(jac.parentSubset == (1:jac.dim(1))'))
            if isa(acc, 'DiagonalJacobian')
                % NB currently not checking subset here - bug
                jac = mexDiscreteDivergenceJac(acc.diagonal, jac.diagonal, N, prelim.facePos, prelim.faces, prelim.cells, prelim.cellIndex);
            else
                jac = acc + mexDiscreteDivergenceJac([], jac.diagonal, N, prelim.facePos, prelim.faces, prelim.cells, prelim.cellIndex);
            end
        else
            jac = acc + sortIx.C*jac.sparse();
        end
    end
end
