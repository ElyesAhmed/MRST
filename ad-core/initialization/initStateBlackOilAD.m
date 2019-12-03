function [state, pressures] = initStateBlackOilAD(model, regions, varargin)
%Undocumented Utility Function

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

    opt = struct('pressure', []);
    opt = merge_options(opt, varargin{:});
    
    if ~iscell(regions)
        regions = {regions};
    end
    
    [rs, rv] = deal(0);
    vapoil = isprop(model, 'vapoil') && model.vapoil;
    disgas = isprop(model, 'disgas') && model.disgas;
    compositional = isprop(model, 'EOSModel');

    G = model.G;
    nph = sum(model.getActivePhases());
    state = struct('pressure', zeros(G.cells.num, 1), 's', zeros(G.cells.num, nph));
    if disgas
        state.rs = zeros(G.cells.num, 1);
    end
    if vapoil
        state.rv = zeros(G.cells.num, 1);
    end
    if compositional
        ncomp = model.EOSModel.fluid.getNumberOfComponents();
        state.components = zeros(G.cells.num, ncomp);
        state.T = zeros(G.cells.num, 1);
    end
    watIx = model.getPhaseIndex('W');
    oilIx = model.getPhaseIndex('O');
    gasIx = model.getPhaseIndex('G');

    pressures = zeros(G.cells.num, nph);
    touched = false(G.cells.num, 1);
    for regNo = 1:numel(regions)
        region = regions{regNo};
        if ischar(region.cells)
            assert(numel(regions) == 1)
            region.cells = (1:model.G.cells.num)';
        end
        cells = region.cells;
        assert(~any(touched(cells)), 'Multiple regions defined in same cells.');
        touched(cells) = true;
        
        if isempty(opt.pressure)
            p = initializeEquilibriumPressures(model, region);
        else
            p = opt.pressure(cells, :);
        end
        z = model.G.cells.centroids(cells, 3);
        
        s = initializeEquilibriumSaturations(model, region, p);
        state.s(cells, :) = s;
        pressures(cells, :) = p;
        
        % Evalaute rel. perm.
        sat = cell(1, nph);
        pc = cell(1, nph);
        for i = 1:nph
            sat{i} = state.s(cells, i);
            pc{i} = region.pc_sign(i)*region.pc_functions{i}(state.s(cells, i));
        end
        kr = s;
        numberOfMobile = sum(kr > 0, 2);
        maxSat = max(s, [], 2);
        singlePhaseMobile = numberOfMobile <= 1;
        
        toOil = true(size(p, 1), 1);
        if model.gas
            % If only gas is mobile, set oil pressure to the gas hydrostatic 
            % pressure minus the capillary pressure
            onlyGas = (kr(:, gasIx) > 0 & singlePhaseMobile) |...
                       (s(:, gasIx) == maxSat & numberOfMobile == 0);

            toOil(onlyGas) = false;
            state.pressure(cells(onlyGas)) = p(onlyGas, gasIx) - pc{gasIx}(onlyGas);
            if disgas
                po = p(:, oilIx);
                if iscell(model.fluid.rsSat)
                    rsSatF = model.fluid.rsSat{region.pvt_region};
                else
                    rsSatF = model.fluid.rsSat;
                end
                rsMax = rsSatF(po);
                rs = region.rs(po, z);
                rs(rs > rsMax) = rsMax(rs > rsMax);
                sg = s(:, gasIx);
                rs(sg > 0) = rsMax(sg > 0);
                state.rs(cells) = rs;             
            end
        end
        if model.oil
            if vapoil
                pg = p(:, gasIx);
                if iscell(model.fluid.rvSat)
                    rvSatF = model.fluid.rvSat{region.pvt_region};
                else
                    rvSatF = model.fluid.rvSat;
                end
                rvMax = rvSatF(po);
                rv = region.rv(pg, z);
                rv(rv > rvMax) = rvMax(rv > rvMax);
                so = s(:, oilIx);
                rv(so > 0) = rvMax(so > 0);
                state.rv(cells) = rv;
            end
        end
        if model.water
            % onlyWat = state.s(cells, watIx) == 1;
            onlyWat = (kr(:, watIx) > 0 & singlePhaseMobile) | ...
                      (s(:, watIx) == maxSat & numberOfMobile == 0);
            toOil(onlyWat) = false;
            state.pressure(cells(onlyWat)) = p(onlyWat, watIx) - pc{watIx}(onlyWat);
        end
        if model.oil
            state.pressure(cells(toOil)) = p(toOil, oilIx);
        end
        
        if compositional
            if isfield(region, 'z')
                state.components(cells, :) = region.z(state.pressure(cells), z);
            end
            if isfield(region, 'T')
                state.T(cells, :) = region.T(state.pressure(cells), z);
            end
        end
    end
    if ~all(touched)
        warning('Regions did not cover all cells. Model only partially initialized.');
    end
end
