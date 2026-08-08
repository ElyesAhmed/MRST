function state = runH2StoragePhreeqcTimestepCoupling(model, state, varargin)
% Equilibrate every H2Storage1D cell after a converged MRST timestep.
%
% The coupling is intentionally operator-split. It maps PHREEQC aqueous
% carbonate/sulfate/Ca/Mg, pH, and finite mineral inventories into state
% for the next timestep. When the explicit
% phreeqcConservativeCarbonTransfer option is enabled, the change in
% PHREEQC total inorganic carbon is also mapped to the EOS CO2 inventory
% and the EOS state is flashed again. Fixed brine ions remain inputs only;
% this is not a fully coupled geochemical transport model.

% Accept the current timestep for a uniform post-convergence dispatch
% interface. The standard equilibrium backend deliberately does not use it.
assert(numel(varargin) <= 1, ...
    'runH2StoragePhreeqcTimestepCoupling accepts at most one optional dt argument.');
if ~isempty(varargin)
    validateattributes(varargin{1}, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'dt');
end

assert(model.phreeqcTimestepCoupling, ...
    'PHREEQC timestep coupling is not enabled on this model.');
assert(strcmp(model.phreeqcBackend, 'standard-pitzer'), ...
    ['runH2StoragePhreeqcTimestepCoupling only supports ', ...
     'phreeqcBackend=''standard-pitzer''.']);
assert(model.carbonateBuffer && model.sulfateReduction, ...
    ['PHREEQC timestep coupling requires the HCO3 and SO4 aqueous ', ...
     'tracers (carbonateBuffer and sulfateReduction).']);
assert(exist('IPhreeqc', 'file') == 2, ...
    ['PHREEQC timestep coupling requires PhreeqcMatlab. Run ', ...
     'startup_Phreeqc before simulating.']);

opt = getCouplingOptions(model);
databaseFile = opt.databaseFile;
if isempty(databaseFile)
    assert(exist('database_file', 'file') == 2, ...
        ['Provide phreeqcDatabaseFile or run startup_Phreeqc before ', ...
         'simulating with phreeqcTimestepCoupling=true.']);
    databaseFile = database_file('pitzer.dat');
end
assert(isfile(databaseFile), 'PHREEQC database not found: %s', databaseFile);
databaseName = lower(string(databaseFile));
assert(contains(databaseName, 'pitzer'), ...
    ['PHREEQC timestep coupling requires a standard PHREEQC Pitzer ', ...
     'database (normally pitzer.dat).']);
assert(~contains(databaseName, 'pitzer_modified') && ~contains(databaseName, 'ugfact'), ...
    ['PHREEQC timestep coupling only supports standard PHREEQC Pitzer ', ...
     'databases; UGFACT PITZER_Modified.DAT is incompatible.']);

required = {'tracerHCO3', 'tracerSO4', 'tracerCa', 'tracerMg', 'phreeqcPH', ...
    'phreeqcMineralDolomite', 'phreeqcMineralCalcite', ...
    'phreeqcMineralBrucite', 'phreeqcMineralQuartz'};
for i = 1:numel(required)
    assert(isfield(state, required{i}), ...
        'PHREEQC timestep state is missing %s. Initialize it before simulation.', required{i});
end

nc = model.G.cells.num;
rhoWater = opt.waterDensity;
waterMass = getWaterMass(model, state, rhoWater);
assert(all(waterMass > 0), ...
    'PHREEQC timestep coupling requires positive liquid water mass in every cell.');

% tracerHCO3 is the transported non-CO2 DIC inventory, not necessarily
% the HCO3- species. This retains carbonate species PHREEQC keeps outside
% the EOS while the actual HCO3- molality is stored separately below.
nonvolatileCarbon = max(value(state.tracerHCO3)./rhoWater, 0);
so4 = max(value(state.tracerSO4)./rhoWater, 0);
[ca, mg] = deal(max(value(state.tracerCa)./rhoWater, 0), ...
    max(value(state.tracerMg)./rhoWater, 0));
[co2, xWater] = getAqueousCarbon(model, state);
co2 = max(co2.*55.508./max(xWater, 1e-12), 0);
totalCarbon = asCellVector(nonvolatileCarbon + co2, nc, 'aqueous carbon');
so4 = asCellVector(so4, nc, 'sulfate');
ca = asCellVector(ca, nc, 'calcium');
mg = asCellVector(mg, nc, 'magnesium');
pH = asCellVector(value(state.phreeqcPH), nc, 'pH');
temperature = asCellVector(value(state.T) - 273.15, nc, 'temperature');
% EQUILIBRIUM_PHASES amounts are moles, as are the mineral state fields.
% Unlike dissolved species, they must not be converted to molality.
dolomite = asCellVector(max(value(state.phreeqcMineralDolomite), 0), ...
    nc, 'Dolomite inventory');
calcite = asCellVector(max(value(state.phreeqcMineralCalcite), 0), ...
    nc, 'Calcite inventory');
brucite = asCellVector(max(value(state.phreeqcMineralBrucite), 0), ...
    nc, 'Brucite inventory');
quartz = asCellVector(max(value(state.phreeqcMineralQuartz), 0), ...
    nc, 'Quartz inventory');
state.phreeqcInputTotalCarbon = totalCarbon;
state.phreeqcInputMineralDolomite = dolomite;

input = buildInput(opt, nc, temperature, pH, totalCarbon, so4, ca, mg, ...
    dolomite, calcite, brucite, quartz, waterMass);
iph = IPhreeqc();
iph = iph.CreateIPhreeqc();
cleanup = onCleanup(@() iph.DestroyIPhreeqc()); %#ok<NASGU>

status = iph.LoadDatabase(databaseFile);
assert(status == 0, 'PHREEQC could not load database: %s', databaseFile);
iph.SetOutputStringOn(1);
iph.SetSelectedOutputStringOn(1);
status = iph.RunString(input);
assert(status == 0, 'PHREEQC timestep equilibrium failed:\n%s', iph.GetErrorString());

[headers, values] = parseSelectedOutput(iph.GetSelectedOutputString());
assert(size(values, 1) >= 2*nc, ...
    'PHREEQC returned %d selected-output rows for %d cells.', size(values, 1), nc);
% Each cell produces an initial-solution row followed by its equilibrated
% USE solution/EQUILIBRIUM_PHASES row. Retain only the latter; selecting
% the final nc rows mixes both types whenever multiple cells are present.
values = values(end - 2*nc + 2:2:end, :);
assert(size(values, 1) == nc, ...
    'PHREEQC selected-output final-row extraction failed for %d cells.', nc);

finalPH = getSelectedValue(headers, values, 'pH');
finalCarbon = max(getSelectedValue(headers, values, 'C_4_'), 0);
finalSulfate = max(getSelectedValue(headers, values, 'S_6_'), 0);
finalCa = max(getSelectedValue(headers, values, 'Ca'), 0);
finalMg = max(getSelectedValue(headers, values, 'Mg'), 0);
finalHCO3 = max(getSelectedValue(headers, values, 'm_HCO3_'), 0);
finalCO2 = max(getSelectedValue(headers, values, 'm_CO2'), 0);
finalDolomite = max(getSelectedValue(headers, values, ...
    matlab.lang.makeValidName(opt.dolomitePhase)), 0);
finalCalcite = max(getSelectedValue(headers, values, ...
    matlab.lang.makeValidName(opt.calcitePhase)), 0);
finalBrucite = max(getSelectedValue(headers, values, ...
    matlab.lang.makeValidName(opt.brucitePhase)), 0);
finalQuartz = max(getSelectedValue(headers, values, ...
    matlab.lang.makeValidName(opt.quartzPhase)), 0);
finalDolomiteSI = getSelectedValue(headers, values, ...
    matlab.lang.makeValidName(['si_', opt.dolomitePhase]));
finalCalciteSI = getSelectedValue(headers, values, ...
    matlab.lang.makeValidName(['si_', opt.calcitePhase]));
finalBruciteSI = getSelectedValue(headers, values, ...
    matlab.lang.makeValidName(['si_', opt.brucitePhase]));

state.tracerHCO3 = max(finalCarbon - finalCO2, 0).*rhoWater;
state.tracerSO4 = finalSulfate.*rhoWater;
state.tracerCa = finalCa.*rhoWater;
state.tracerMg = finalMg.*rhoWater;
state.phreeqcPH = finalPH;
state.phreeqcHCO3Molality = finalHCO3;
carbonatePka1 = finalPH - log10(finalHCO3./max(finalCO2, 1e-30));
if isfield(state, 'phreeqcCarbonatePka1')
    previousPka1 = asCellVector(value(state.phreeqcCarbonatePka1), nc, ...
        'previous carbonate pKa1');
else
    previousPka1 = repmat(model.carbonateBufferPka1, nc, 1);
end
carbonatePka1(~isfinite(carbonatePka1)) = previousPka1(~isfinite(carbonatePka1));
state.phreeqcCarbonatePka1 = carbonatePka1;
state.phreeqcCO2Molality = finalCO2;
state.phreeqcTotalCarbon = finalCarbon;
state.phreeqcMineralDolomite = finalDolomite;
state.phreeqcMineralCalcite = finalCalcite;
state.phreeqcMineralBrucite = finalBrucite;
state.phreeqcMineralQuartz = finalQuartz;
state.phreeqcDolomiteSI = finalDolomiteSI;
state.phreeqcCalciteSI = finalCalciteSI;
state.phreeqcBruciteSI = finalBruciteSI;

if model.phreeqcConservativeCarbonTransfer
    initialCarbonMoles = totalCarbon.*waterMass;
    finalCarbonMoles = finalCarbon.*waterMass;
    initialNonvolatileCarbonMoles = nonvolatileCarbon.*waterMass;
    finalNonvolatileCarbonMoles = max(finalCarbon - finalCO2, 0).*waterMass;

    % Non-CO2 DIC is represented by its own MRST tracer. Its PHREEQC
    % update already changes that carbon inventory, so transfer only the
    % remaining DIC change to the EOS CO2 component.
    carbonTransferMoles = (finalCarbonMoles - initialCarbonMoles) - ...
        (finalNonvolatileCarbonMoles - initialNonvolatileCarbonMoles);
    assert(all(isfinite(carbonTransferMoles)), ...
        'PHREEQC total-inorganic-carbon transfer contains non-finite values.');

    state = transferCarbonToEOSCO2(model, state, carbonTransferMoles);
    state.phreeqcCarbonTransferMoles = carbonTransferMoles;
    if isfield(state, 'phreeqcCumulativeCarbonTransferMoles')
        cumulative = asCellVector(value(state.phreeqcCumulativeCarbonTransferMoles), ...
            nc, 'cumulative PHREEQC carbon transfer');
    else
        cumulative = zeros(nc, 1);
    end
    state.phreeqcCumulativeCarbonTransferMoles = cumulative + carbonTransferMoles;
else
    % pH, DIC, Ca, and Mg changed outside the Newton solve. Drop cached
    % state functions so the next timestep evaluates the transferred state.
    state = clearStateFunctionCaches(model, state);
end
end

function state = transferCarbonToEOSCO2(model, state, carbonTransferMoles)
% Add PHREEQC's net mineral-carbon change to the total EOS CO2 inventory.
%
% State.components holds overall mole fractions, not component inventories.
% Recover each cell's EOS moles from pore volume, phase saturations, and
% phase molar densities before adding the one-carbon-per-CO2 transfer.

names = model.EOSModel.CompositionalMixture.names;
iCO2 = find(strcmp(names, 'CO2'), 1);
assert(~isempty(iCO2), ...
    'Conservative PHREEQC carbon transfer requires an EOS CO2 component.');

componentMoles = getEOSComponentMoles(model, state);
updatedComponentMoles = componentMoles;
updatedComponentMoles(:, iCO2) = updatedComponentMoles(:, iCO2) + carbonTransferMoles;
assert(all(updatedComponentMoles(:, iCO2) >= 0), ...
    ['PHREEQC carbon transfer would make the EOS CO2 inventory negative. ', ...
     'The PHREEQC carbon decrease exceeds available EOS CO2.']);

updatedTotalMoles = sum(updatedComponentMoles, 2);
assert(all(isfinite(updatedTotalMoles) & updatedTotalMoles > 0), ...
    'EOS total moles are invalid after the PHREEQC carbon transfer.');
state.components = bsxfun(@rdivide, updatedComponentMoles, updatedTotalMoles);
assert(all(isfinite(state.components(:))) && ...
    all(abs(sum(state.components, 2) - 1) < 100*eps), ...
    'Conservative PHREEQC carbon transfer produced invalid EOS compositions.');

% A changed overall composition invalidates cached PVT/flow properties.
% Reflash before the next timestep and retain fresh requested outputs.
state = clearStateFunctionCaches(model, state);
model = updateEOSSalinityForReflash(model, state);
state = model.computeFlash(state, inf);
state = refreshOutputStateFunctions(model, state);
end

function componentMoles = getEOSComponentMoles(model, state)
nc = model.G.cells.num;
poreVolume = asCellVector(value(model.PVTPropertyFunctions.get( ...
    model, state, 'PoreVolume')), nc, 'pore volume');
s = value(state.s);
liquid = model.getLiquidIndex();
vapor = model.getVaporIndex();
if iscell(s)
    sL = s{liquid};
    sV = s{vapor};
else
    sL = s(:, liquid);
    sV = s(:, vapor);
end
sL = asCellVector(sL, nc, 'liquid saturation');
sV = asCellVector(sV, nc, 'vapor saturation');

propmodel = model.EOSModel.PropertyModel;
rhoL = propmodel.computeMolarDensity(model.EOSModel, value(state.pressure), ...
    value(state.x), value(state.Z_L), value(state.T), true);
rhoV = propmodel.computeMolarDensity(model.EOSModel, value(state.pressure), ...
    value(state.y), value(state.Z_V), value(state.T), false);
rhoL = asCellVector(value(rhoL), nc, 'liquid molar density');
rhoV = asCellVector(value(rhoV), nc, 'vapor molar density');
totalMoles = poreVolume.*(sL.*rhoL + sV.*rhoV);
assert(all(isfinite(totalMoles) & totalMoles > 0), ...
    'Conservative PHREEQC carbon transfer requires positive EOS moles per cell.');

components = value(state.components);
assert(ismatrix(components) && size(components, 1) == nc && ...
    size(components, 2) == model.EOSModel.getNumberOfComponents(), ...
    'EOS components must be a cell-by-component matrix for carbon transfer.');
componentMoles = bsxfun(@times, totalMoles, components);
assert(all(isfinite(componentMoles(:)) & componentMoles(:) >= 0), ...
    'Conservative PHREEQC carbon transfer requires non-negative EOS inventories.');
end

function state = clearStateFunctionCaches(model, state)
groups = model.getStateFunctionGroupings();
for i = 1:numel(groups)
    name = groups{i}.getStateFunctionContainerName();
    if isfield(state, name)
        state = rmfield(state, name);
    end
end
end

function state = refreshOutputStateFunctions(model, state)
outputs = unique(model.OutputStateFunctions, 'stable');
for i = 1:numel(outputs)
    [~, state] = model.getProp(state, outputs{i});
end
end

function model = updateEOSSalinityForReflash(model, state)
if model.sulfateReduction && isa(model.EOSModel, 'SoreideWhitsonEos')
    assert(isfield(state, 'tracerHS') && isfield(state, 'h2sDissolvedLag'), ...
        'PHREEQC carbon transfer requires HS and lagged dissolved-H2S state.');
    model.EOSModel = model.EOSModel.enablesrb_coupling( ...
        value(state.tracerSO4), value(state.tracerHS), ...
        value(state.tracerHS) + value(state.h2sDissolvedLag), value(state.T));
end
end

function input = buildInput(opt, nc, temperature, pH, carbon, sulfate, ca, mg, ...
        dolomite, calcite, brucite, quartz, waterMass)
input = sprintf([ ...
    'SELECTED_OUTPUT 1\n' ...
    '    -reset false\n' ...
    '    -pH true\n' ...
    '    -totals C(4) S(6) Ca Mg\n' ...
    '    -molalities CO2 HCO3-\n' ...
    '    -equilibrium_phases %s %s %s %s\n' ...
    '    -saturation_indices %s %s %s\n' ...
    'END\n'], opt.dolomitePhase, opt.calcitePhase, opt.brucitePhase, opt.quartzPhase, ...
    opt.dolomitePhase, opt.calcitePhase, opt.brucitePhase);
for cellNo = 1:nc
    input = [input, sprintf([ ... %#ok<AGROW>
        'SOLUTION %d\n' ...
        '    units mol/kgw\n' ...
        '    -water %.15g\n' ...
        '    temp %.15g\n' ...
        '    pH %.15g\n' ...
        '    Na %.15g\n' ...
        '    Ca %.15g\n' ...
        '    Mg %.15g\n' ...
        '    Cl %.15g\n' ...
        '    S(6) %.15g as SO4\n' ...
        '    C(4) %.15g as HCO3\n' ...
        '    Si %.15g\n' ...
        'EQUILIBRIUM_PHASES %d\n' ...
        '    %s 0 %.15g\n' ...
        '    %s 0 %.15g\n' ...
        '    %s 0 %.15g\n' ...
        '    %s 0 %.15g\n' ...
        'USE solution %d\n' ...
        'USE equilibrium_phases %d\n' ...
        'END\n'], ...
        cellNo, waterMass(cellNo), temperature(cellNo), pH(cellNo), ...
        optionAt(opt.Na, cellNo), ca(cellNo), mg(cellNo), optionAt(opt.Cl, cellNo), ...
        sulfate(cellNo), carbon(cellNo), optionAt(opt.Si, cellNo), ...
        cellNo, opt.dolomitePhase, dolomite(cellNo), ...
        opt.calcitePhase, calcite(cellNo), ...
        opt.brucitePhase, brucite(cellNo), ...
        opt.quartzPhase, quartz(cellNo), cellNo, cellNo)];
end
end

function [co2, xWater] = getAqueousCarbon(model, state)
names = model.EOSModel.CompositionalMixture.names;
iCO2 = find(strcmp(names, 'CO2'), 1);
iWater = find(strcmp(names, 'H2O'), 1);
assert(~isempty(iCO2) && ~isempty(iWater), ...
    'PHREEQC timestep coupling requires CO2 and H2O EOS components.');
x = value(state.x);
if iscell(x)
    co2 = x{iCO2};
    xWater = x{iWater};
else
    co2 = x(:, iCO2);
    xWater = x(:, iWater);
end
end

function waterMass = getWaterMass(model, state, rhoWater)
poreVolume = value(model.PVTPropertyFunctions.get(model, state, 'PoreVolume'));
s = value(state.s);
liquid = model.getLiquidIndex();
if iscell(s)
    sL = s{liquid};
else
    sL = s(:, liquid);
end
waterMass = poreVolume.*max(sL, 1e-12).*rhoWater;
end

function valueAtCell = optionAt(value, cellNo)
if isscalar(value)
    valueAtCell = value;
else
    valueAtCell = value(cellNo);
end
end

function values = asCellVector(values, nc, name)
if isscalar(values)
    values = repmat(values, nc, 1);
else
    values = values(:);
    assert(numel(values) == nc, ...
        'PHREEQC timestep coupling expected %d %s values, got %d.', ...
        nc, name, numel(values));
end
end

function opt = getCouplingOptions(model)
opt = struct( ...
    'databaseFile', model.phreeqcDatabaseFile, ...
    'temperature', 60, ...
    'waterDensity', 1000, ...
    'Na', 2.865, ...
    'Cl', 3.655, ...
    'Si', 9.723e-5, ...
    'dolomitePhase', 'Dolomite', ...
    'calcitePhase', 'Calcite', ...
    'brucitePhase', 'Brucite', ...
    'quartzPhase', 'Quartz');
configured = model.phreeqcCouplingOptions;
names = fieldnames(opt);
for i = 1:numel(names)
    name = names{i};
    if isfield(configured, name)
        opt.(name) = configured.(name);
    end
end
end

function [headers, values] = parseSelectedOutput(selectedOutput)
lines = splitlines(string(selectedOutput));
lines = strtrim(lines(strlength(strtrim(lines)) > 0));
assert(numel(lines) >= 2, 'PHREEQC selected output contains no data rows.');
headers = regexp(char(lines(1)), '\s+', 'split');
headers = cellfun(@matlab.lang.makeValidName, headers, 'UniformOutput', false);
values = zeros(numel(lines) - 1, numel(headers));
for i = 2:numel(lines)
    row = sscanf(char(lines(i)), '%f').';
    assert(numel(row) == numel(headers), ...
        'Could not parse PHREEQC selected-output row %d.', i - 1);
    values(i - 1, :) = row;
end
end

function selected = getSelectedValue(headers, values, name)
index = find(strcmp(headers, name), 1);
assert(~isempty(index), 'PHREEQC selected output is missing column %s.', name);
selected = values(:, index);
end
