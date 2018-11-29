function v = faceAverage(N, v, useMex)
    % Face average operator for the NewAD library
    if isa(v, 'NewAD')
        v.val = 0.5*sum(v.val(N), 2);
        v.jac = cellfun(@(x) avgJac(x, N, useMex), v.jac, 'UniformOutput', false);
    else
        v = 0.5*(v(N(:, 1), :) + v(N(:, 2), :));
    end
end

function jac = avgJac(jac, N, useMex)
    if issparse(jac)
        if any(jac(:))
            jac = 0.5*(jac(N(:, 1), :) + jac(N(:, 2), :));
        else
            jac = sparse([], [], [], size(N, 1), size(jac, 2));
        end
    elseif jac.isZero
        jac = jac.toZero(size(N, 1));
    else
        if useMex
            diagonal = mexFaceAverageDiagonalJac(jac.diagonal, N);
        else
            diagonal = 0.5*jac.diagonal(N, :);
        end
        if isempty(jac.subset)
            map = N;
        else
            map = jac.subset(N);
        end
        jac = DiagonalSubset(diagonal, jac.dim, map);
    end

end