function [H2_rate_per_cell, H2_cum_per_cell] = computeH2Consumption(states, schedule, model, idxReaction)

bf = model.biochemFluid;
gam_H2 = bf.gamrH2(idxReaction);
psi_max = bf.Psigrowthmax(idxReaction);
l_H2 = bf.alphaH2(idxReaction);
l_sub = bf.alphasub(idxReaction);
Y_H2 = bf.Y_H2(idxReaction);          % cells/mol H2 or mol biomass/mol H2
n0 = bf.nbactMax(idxReaction);        % cells/m3
L_ix = model.getLiquidIndex();
% --- Determine if Y_H2 is in cells/mol or mol biomass/mol ---
% If Y_H2 is large (>1e10), it's in cells/mol; convert to mol biomass/mol
if Y_H2 > 1e10
    cells_per_mol_biomass = 3.333e12;     % typical conversion factor
    Y_H2 = Y_H2 / cells_per_mol_biomass;  % now in mol biomass/mol H2
    % Also convert n0 from cells/m3 to mol biomass/m3
    n0 = n0 / cells_per_mol_biomass;  % now in mol biomass/m3
end

isSRB = strcmp(bf.metabolicReaction(idxReaction), 'SulfateReducingBacteria');

eosNames = model.EOSModel.CompositionalMixture.names;
idx_H2 = find(strcmp(eosNames, 'H2'));
mc_H2 = model.EOSModel.CompositionalMixture.molarMass(idx_H2);
if isSRB
    idx_sub = [];
else
    idx_sub = find(strcmp(eosNames, 'CO2'));
end

volumes = model.G.cells.volumes;
nCells = numel(volumes);
nSteps = numel(states);
H2_rate_per_cell = zeros(nCells, nSteps);
H2_cum_per_cell = zeros(nCells, nSteps);

for t = 1:nSteps
    state = states{t};

    if t == 1
        dt = schedule.step.val(1);
    else
        dt = schedule.step.val(t);
    end

    rhoL = state.FlowProps.ComponentPhaseDensity{L_ix};
    x_H2 = getNumeric(state, 'x', idx_H2);

    if isSRB
        SO4_conc = state.tracerSO4;
        x_sub = SO4_conc ./ rhoL;
    else
        x_sub = getNumeric(state, 'x', idx_sub);
    end

    n_tilde = getNumeric(state, 'nbact', idxReaction);
    s_l = getNumeric(state, 's', L_ix);

    % --- Porosity handling (supports bioclogging) ---
    if isfield(state, 'porosity')
        Phi = state.porosity;
    else
        % If porosity is not stored in state, evaluate the function handle
        if isa(model.rock.poro, 'function_handle')
            p = state.pressure;
            n_tilde = getNumeric(state, 'nbact', idxReaction);
            nbactArray = model.extractBactValues(state.nbact);
            Phi = model.rock.poro(p, nbactArray{:});
        else
            Phi = model.rock.poro;
        end
    end

    % Monod terms and growth rate
    monod_H2 = x_H2 ./ (l_H2 + x_H2);
    monod_sub = x_sub ./ (l_sub + x_sub);
    psi_growth = psi_max .* monod_H2 .* monod_sub;

    % nbact is the raw, un-rescaled state variable (see BacterialMass.m /
    % BactConvertionRate.m: bmass = pv*sL*rhoL*nbact, no extra nbactMax
    % applied to nbact itself -- nbactMax only enters once, below, exactly
    % as in BactConvertionRate's gamma_norm).
    n_actual = n_tilde;

    % Source term (mol/m3/s). Matches BiochemistryModel.m's actual use of
    % this quantity: BactConvertionRate returns qbiot_H2 [kg/s] with an
    % uncancelled rhoL baked in from bmass; BiochemistryModel.m (line
    % ~301) divides by rhoL when inserting the source into the component
    % equation -- that is where the density actually cancels, not inside
    % BactConvertionRate.m. Converting that mass rate to mol/s divides by
    % mc_H2, which cancels the mc_H2 introduced via the stoichiometric
    % normalization (nbactMax*gam_H2*mc_H2/abs(gam_H2)), so the sign
    % flip -gam_H2/abs(gam_H2) reduces to +-1 and mc_H2/rhoL both drop out.
    r_H2 = Phi .* psi_growth .* n_actual .* s_l .* n0 ./ Y_H2;

    % Rate per cell (mol/s)
    H2_rate_per_cell(:, t) = r_H2 .* volumes;

    % Cumulative per cell (mol)
    if t == 1
        H2_cum_per_cell(:, t) = H2_rate_per_cell(:, t) * dt;
    else
        H2_cum_per_cell(:, t) = H2_cum_per_cell(:, t-1) + H2_rate_per_cell(:, t) * dt;
    end
end
end

function vec = getNumeric(state, fld, idx)
data = state.(fld);
if iscell(data)
    vec = data{idx};
else
    vec = data(:, idx);
end
end
