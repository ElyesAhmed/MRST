function solver = setUpGeothermalLinearSolverAD_WB(model, varargin)

    opt = struct('useCPR', true);
    opt = merge_options(opt, varargin{:});

    if ~opt.useCPR
        solver = AMGCLSolverAD();
        solver.setPreconditioner('relaxation');
        solver.setRelaxation('ilu0');
        solver.setSolver('gmres');
    else
        solver = AMGCL_CPRSolverAD();
        solver.decoupling = 'quasiimpes';
    end
    solver.maxIterations = 200;
    solver.tolerance = 1e-4;
    
    % Set block size
    bz = 1 + model.submodels.Reservoir.thermal;
    solver.amgcl_setup.block_size = bz;
    % Set equation/variable order for reservor dofs
    ncr = model.submodels.Reservoir.G.cells.num;
    rorder = repmat(1:ncr, bz, 1) + (0:bz-1)'.*ncr;
    rorder = rorder(:);
    % Set equation/variable order for well dofs
    ncw = model.submodels.Wellbore.G.cells.num;
    worder = repmat(1:ncw, bz, 1) + (0:bz-1)'.*ncw + bz*ncr;
    worder = worder(:);
    % Set equation/variable order
    order = [rorder; worder];
    [solver.equationOrdering, solver.variableOrdering] = deal(order);
    % Set cells to keep in Schur complement
    cells = model.submodels.Wellbore.G.cells.type == 0;
    if isfield(model.submodels.Wellbore.G.cells, 'hybrid')
        cells = cells & model.submodels.Wellbore.G.cells.hybrid == 0;
    end
    ncw = nnz(cells);
    solver.keepNumber = bz*(ncr + ncw);
    solver.amgcl_setup.cell_size = ncr + ncw;
    solver.reduceToCell = false;
    
end