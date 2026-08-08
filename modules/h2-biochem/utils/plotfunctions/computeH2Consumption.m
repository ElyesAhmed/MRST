function [H2_rate_per_cell, H2_cum_per_cell] = computeH2Consumption(states, schedule, model, idxReaction)
% Compute reaction-specific H2 consumption from the model source terms.

bf = model.biochemFluid;
nReactions = bf.nbioreact;
assert(idxReaction >= 1 && idxReaction <= nReactions, ...
    'idxReaction must identify an existing biochemical reaction.');
assert(numel(states) == numel(schedule.step.val), ...
    'states and schedule must contain the same number of timesteps.');

nCells = model.G.cells.num;
nSteps = numel(states);
H2_rate_per_cell = zeros(nCells, nSteps);
H2_cum_per_cell = zeros(nCells, nSteps);

for t = 1:nSteps
    state = states{t};

    if model.isMrstMonodComPhreeqcBackend() && ...
            isfield(state, 'phreeqcMrstMonodComCumulativeH2ConsumptionMoles')
        cumulative = state.phreeqcMrstMonodComCumulativeH2ConsumptionMoles(:, idxReaction);
        if t == 1
            increment = cumulative;
        else
            previous = states{t - 1}.phreeqcMrstMonodComCumulativeH2ConsumptionMoles(:, idxReaction);
            increment = cumulative - previous;
        end
        H2_rate_per_cell(:, t) = increment./schedule.step.val(t);
        H2_cum_per_cell(:, t) = cumulative;
        continue;
    end

    assert(isfield(state, 'h2ConsumptionRate'), ...
        ['Saved states do not contain the converged H2 source rate. ', ...
         'Rerun the simulation before post-processing H2 consumption.']);
    H2_rate_per_cell(:, t) = state.h2ConsumptionRate(:, idxReaction);

    if t == 1
        H2_cum_per_cell(:, t) = H2_rate_per_cell(:, t).*schedule.step.val(t);
    else
        H2_cum_per_cell(:, t) = H2_cum_per_cell(:, t - 1) + ...
            H2_rate_per_cell(:, t).*schedule.step.val(t);
    end
end
end
