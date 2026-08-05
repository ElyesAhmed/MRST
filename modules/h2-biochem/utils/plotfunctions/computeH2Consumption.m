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

    % This is the same evaluated kinetic rate and bacterial mass used by
    % BactConvertionRate. BacterialMass already contains pore volume,
    % saturation, density, and any dynamic porosity from bio-clogging.
    psiGrowth = model.FlowPropertyFunctions.get(model, state, 'PsiGrowthRate');
    bacterialMass = model.PVTPropertyFunctions.get(model, state, 'BacterialMass');
    density = model.PVTPropertyFunctions.get(model, state, 'Density');
    psiGrowth = reactionValue(psiGrowth, idxReaction);
    bacterialMass = reactionValue(bacterialMass, idxReaction);
    rhoL = phaseValue(density, model.getLiquidIndex());

    % BiochemistryModel inserts BactConvRate divided by the liquid density
    % into the component equations. Apply the same scaling here so this
    % post-processing reports the H2 consumption represented by the solve.
    H2_rate_per_cell(:, t) = bf.nbactMax(idxReaction).*psiGrowth.* ...
        bacterialMass./(bf.Y_H2(idxReaction).*rhoL);

    if t == 1
        H2_cum_per_cell(:, t) = H2_rate_per_cell(:, t).*schedule.step.val(t);
    else
        H2_cum_per_cell(:, t) = H2_cum_per_cell(:, t - 1) + ...
            H2_rate_per_cell(:, t).*schedule.step.val(t);
    end
end
end

function values = reactionValue(data, idxReaction)
if iscell(data)
    values = data{idxReaction};
else
    values = data(:, idxReaction);
end
end

function values = phaseValue(data, idxPhase)
if iscell(data)
    values = data{idxPhase};
else
    values = data(:, idxPhase);
end
end
