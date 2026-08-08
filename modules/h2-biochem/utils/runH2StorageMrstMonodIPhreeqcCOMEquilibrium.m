function state = runH2StorageMrstMonodIPhreeqcCOMEquilibrium(model, state, varargin)
% Equilibrate PHREEQC_Modified.DAT chemistry without PHREEQC kinetics.
%
% This is the chemistry half of the mrst-monod-com Picard scheme. MRST has
% already advanced its nbact-based Monod reactions when this function is
% called. Consequently, this input intentionally contains no RATES or
% KINETICS block: PHREEQC only repartitions the complete post-reaction
% inventories among aqueous species, gases, and equilibrium minerals.

runtimeOpt = merge_options(struct('inputOnly', false), varargin{:});
validateattributes(runtimeOpt.inputOnly, {'logical'}, {'scalar'}, ...
    mfilename, 'inputOnly');
assert(model.phreeqcTimestepCoupling && model.isMrstMonodComPhreeqcBackend(), ...
    'mrst-monod-com equilibrium chemistry is not enabled on this model.');

opt = getCouplingOptions(model);
validateCouplingOptions(opt);
nc = model.G.cells.num;
waterMass = getWaterMass(model, state, opt.waterDensity);
assert(all(isfinite(waterMass) & waterMass > 0), ...
    'mrst-monod-com requires positive liquid water mass in every cell.');

phase = getPhaseData(model, state, waterMass, opt.waterDensity);
mineral = getMineralData(state, waterMass, nc);
state = storePhreeqcInputDiagnostics(state, phase);
if runtimeOpt.inputOnly
    state.phreeqcMrstMonodComPreReactionInputCell1String = ...
        buildPhreeqcInput(opt, phase, mineral, 1);
    return;
end

assert(ispc, ['mrst-monod-com requires Windows and a registered ', ...
    'IPhreeqcCOM server.']);
assert(exist('actxserver', 'file') == 2 || exist('actxserver', 'builtin') == 5, ...
    ['IPhreeqcCOM activation is unavailable. Install/register the Windows ', ...
     'IPhreeqcCOM server identified by phreeqcComProgId.']);
results = repmat(emptyResult(), nc, 1);
for cellNo = 1:nc
    input = buildPhreeqcInput(opt, phase, mineral, cellNo);
    if cellNo == 1
        state.phreeqcMrstMonodComInputCell1String = input;
    end
    raw = runIPhreeqcCOM(input, opt, cellNo);
    results(cellNo) = parseSelectedOutput(raw, cellNo);
end

state = updateStateFromResults(model, state, opt, phase, results, waterMass);
end

function state = storePhreeqcInputDiagnostics(state, phase)
% Retain the post-MRST, pre-equilibrium PHREEQC basis for Picard audits.
state.phreeqcMrstMonodComInputWaterMass = phase.waterMass;
state.phreeqcMrstMonodComInputH2Moles = phase.gasH2;
state.phreeqcMrstMonodComInputCO2Moles = phase.gasCO2;
state.phreeqcMrstMonodComInputCH4Moles = phase.gasCH4;
state.phreeqcMrstMonodComInputH2SMoles = phase.gasH2S;
state.phreeqcMrstMonodComInputGasVolumeLPerKg = phase.gasVolumeLPerKg;
state.phreeqcMrstMonodComInputPH = phase.pH;
state.phreeqcMrstMonodComInputC4 = phase.C4;
state.phreeqcMrstMonodComInputS6 = phase.S6;
state.phreeqcMrstMonodComInputS2 = phase.S2;
state.phreeqcMrstMonodComInputAcetate = phase.acetate;
state.phreeqcMrstMonodComInputCa = phase.Ca;
state.phreeqcMrstMonodComInputMg = phase.Mg;
end

function opt = getCouplingOptions(model)
opt = struct( ...
    'databaseFile', model.phreeqcDatabaseFile, ...
    'comProgId', model.phreeqcComProgId, ...
    'waterDensity', 1000, ...
    'Na', 2.865, ...
    'K', 0, ...
    'Ca', 0.2857, ...
    'Mg', 0.1144, ...
    'Cl', 3.655, ...
    'Si', 9.723e-5, ...
    'Fe3', 0, ...
    'Fe2', 0);
configured = model.phreeqcCouplingOptions;
for name = fieldnames(opt).'
    if isfield(configured, name{1})
        opt.(name{1}) = configured.(name{1});
    end
end
end

function validateCouplingOptions(opt)
assert(ischar(opt.databaseFile) || ...
    (isstring(opt.databaseFile) && isscalar(opt.databaseFile)), ...
    'mrst-monod-com databaseFile must be a character vector or scalar string.');
databaseFile = char(opt.databaseFile);
assert(~isempty(strtrim(databaseFile)) && isAbsolutePath(databaseFile), ...
    ['mrst-monod-com requires an explicit absolute databaseFile path to ', ...
     'PHREEQC_Modified.DAT.']);
assert(isfile(databaseFile), ...
    'mrst-monod-com PHREEQC database not found: %s', databaseFile);
assert(contains(lower(databaseFile), 'phreeqc_modified.dat'), ...
    ['mrst-monod-com requires PHREEQC_Modified.DAT, not a standard ', ...
     'PHREEQC database: %s'], databaseFile);
assert(ischar(opt.comProgId) || (isstring(opt.comProgId) && isscalar(opt.comProgId)), ...
    'mrst-monod-com comProgId must identify a registered IPhreeqcCOM server.');
assert(~isempty(strtrim(char(opt.comProgId))), ...
    'mrst-monod-com comProgId must identify a registered IPhreeqcCOM server.');
validateattributes(opt.waterDensity, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'waterDensity');
for name = {'Na', 'K', 'Ca', 'Mg', 'Cl', 'Si', 'Fe3', 'Fe2'}
    values = opt.(name{1});
    validateattributes(values, {'numeric'}, {'real', 'finite', 'nonnegative'}, ...
        mfilename, name{1});
end
end

function phase = getPhaseData(model, state, waterMass, rhoWater)
nc = model.G.cells.num;
phase = struct();
phase.waterMass = waterMass;
phase.temperature = asCellVector(value(state.T), nc, 'temperature');
phase.pressureAtm = asCellVector(value(state.pressure)./101325, nc, 'pressure');
phase.pH = getStateVector(state, 'phreeqcPH', nc, model.carbonateBufferPH);
phase.pe = getStateVector(state, 'phreeqcPE', nc, 4);
phase.C4 = max(getStateVector(state, 'tracerHCO3', nc, 0)./rhoWater, 0);
phase.S6 = max(getStateVector(state, 'tracerSO4', nc, 0)./rhoWater, 0);
phase.S2 = max(getStateVector(state, 'tracerHS', nc, 0)./rhoWater, 0);
phase.Ca = max(getStateVector(state, 'tracerCa', nc, 0)./rhoWater, 0);
phase.Mg = max(getStateVector(state, 'tracerMg', nc, 0)./rhoWater, 0);
phase.Na = optionVector(model.phreeqcCouplingOptions, 'Na', 2.865, nc);
phase.K = optionVector(model.phreeqcCouplingOptions, 'K', 0, nc);
phase.Cl = optionVector(model.phreeqcCouplingOptions, 'Cl', 3.655, nc);
phase.Si = optionVector(model.phreeqcCouplingOptions, 'Si', 9.723e-5, nc);
phase.Fe3 = optionVector(model.phreeqcCouplingOptions, 'Fe3', 0, nc);
phase.Fe2 = optionVector(model.phreeqcCouplingOptions, 'Fe2', 0, nc);

phase.componentMoles = getEOSComponentMoles(model, state);
names = model.EOSModel.CompositionalMixture.names;
validateSupportedComponents(names);
phase.gasH2 = getComponentMoles(phase.componentMoles, names, {'H2', 'Hydrogen'});
phase.gasCO2 = getComponentMoles(phase.componentMoles, names, {'CO2', 'CarbonDioxide'});
phase.gasCH4 = getComponentMoles(phase.componentMoles, names, {'C1', 'CH4', 'Methane'});
phase.gasH2S = getComponentMoles(phase.componentMoles, names, {'H2S', 'HydrogenSulfide'});
phase.gasN2 = getComponentMoles(phase.componentMoles, names, {'N2', 'Nitrogen'}, false);
phase.acetate = getComponentMoles(phase.componentMoles, names, ...
    {'CH3COOH', 'AceticAcid', 'Acetate'}, false)./waterMass;

% Supply PHREEQC with complete volatile inventories. The analytical C(4),
% S(6), S(-2), Ca, and Mg totals are separate MRST aqueous tracers, so no
% dissolved EOS component is added to them here.
[~, ~, ~, sV] = getPhaseProperties(model, state);
poreVolume = asCellVector(value(model.PVTPropertyFunctions.get( ...
    model, state, 'PoreVolume')), nc, 'pore volume');
phase.gasVolumeLPerKg = max(poreVolume.*sV.*1000./waterMass, 1e-12);
end

function mineral = getMineralData(state, waterMass, nc)
mineral = struct();
mineral.calcite = getStateVector(state, 'phreeqcMineralCalcite', nc, 0)./waterMass;
mineral.dolomite = getStateVector(state, 'phreeqcMineralDolomite', nc, 0)./waterMass;
mineral.anhydrite = getStateVector(state, 'phreeqcMineralAnhydrite', nc, 0)./waterMass;
mineral.quartz = getStateVector(state, 'phreeqcMineralQuartz', nc, 0)./waterMass;
mineral.goethite = getStateVector(state, 'phreeqcMineralGoethite', nc, 0)./waterMass;
mineral.brucite = getStateVector(state, 'phreeqcMineralBrucite', nc, 0)./waterMass;
mineral.portlandite = getStateVector(state, 'phreeqcMineralPortlandite', nc, 0)./waterMass;
mineral.pyrite = getStateVector(state, 'phreeqcMineralPyrite', nc, 0)./waterMass;
mineral.gypsum = getStateVector(state, 'phreeqcMineralGypsum', nc, 0)./waterMass;
end

function input = buildPhreeqcInput(opt, phase, mineral, cellNo)
totalGas = phase.gasH2(cellNo) + phase.gasCO2(cellNo) + ...
    phase.gasCH4(cellNo) + phase.gasH2S(cellNo) + phase.gasN2(cellNo);
partialPressure = @(moles) phase.pressureAtm(cellNo).*moles./max(totalGas, 1e-30);

input = sprintf([ ...
    'KNOBS\n' ...
    '-iterations 800\n' ...
    '-step_size 30\n' ...
    '-convergence_tolerance 1e-10\n' ...
    'END\n' ...
    'SOLUTION 1\n' ...
    '-pressure %.15g\n' ...
    '-temp %.15g\n' ...
    'pH %.15g #charge\n' ...
    'pe %.15g\n' ...
    'units mol/kgw\n' ...
    'K %.15g\n' ...
    'Na %.15g\n' ...
    'Mg %.15g\n' ...
    'Ca %.15g\n' ...
    'Cl %.15g\n' ...
    'Carbonate(4) %.15g\n' ...
    'Sulfate(6) %.15g\n' ...
    'Sulfide(-2) %.15g\n' ...
    'Fe_tri %.15g\n' ...
    'Fe_di %.15g\n' ...
    'Si %.15g\n' ...
    'Acetate %.15g\n' ...
    '-water 1\n' ...
    'END\n' ...
    'GAS_PHASE 1\n' ...
    '-pressure %.15g\n' ...
    '-temp %.15g\n' ...
    '-fixed_volume\n' ...
    '-volume %.15g\n' ...
    'H2(g) %.15g\n' ...
    'redoxCarbonateO2(g) %.15g\n' ...
    'redoxCH4(g) %.15g\n' ...
    'redoxH2S(g) %.15g\n' ...
    'N2(g) %.15g\n' ...
    'GAS_PHASE_MODIFY 1\n' ...
    '-type 1\n' ...
    '-total_p %.15g\n' ...
    '-volume %.15g\n' ...
    '-component H2(g)\n' ...
    '-moles %.15g\n' ...
    '-component redoxCH4(g)\n' ...
    '-moles %.15g\n' ...
    '-component redoxCarbonateO2(g)\n' ...
    '-moles %.15g\n' ...
    '-component redoxH2S(g)\n' ...
    '-moles %.15g\n' ...
    '-component N2(g)\n' ...
    '-moles %.15g\n' ...
    'END\n' ...
    'EQUILIBRIUM_PHASES 1\n' ...
    'redoxCalcite 0 %.15g\n' ...
    'redoxDolomite 0 %.15g\n' ...
    'redoxAnhydrite 0 %.15g\n' ...
    'Quartz 0 %.15g\n' ...
    'redoxGoethite 0 %.15g\n' ...
    'redoxPyrite 0 %.15g\n' ...
    'Brucite 0 %.15g\n' ...
    'Portlandite 0 %.15g\n' ...
    'redoxGypsum 0 %.15g\n' ...
    'END\n' ...
    'USE solution 1\n' ...
    'USE equilibrium_phases 1\n' ...
    'USE gas_phase 1\n' ...
    '%s'], ...
    phase.pressureAtm(cellNo), phase.temperature(cellNo) - 273.15, ...
    phase.pH(cellNo), phase.pe(cellNo), phase.K(cellNo), phase.Na(cellNo), ...
    phase.Mg(cellNo), phase.Ca(cellNo), phase.Cl(cellNo), phase.C4(cellNo), ...
    phase.S6(cellNo), phase.S2(cellNo), phase.Fe3(cellNo), phase.Fe2(cellNo), ...
    phase.Si(cellNo), phase.acetate(cellNo), phase.pressureAtm(cellNo), ...
    phase.temperature(cellNo) - 273.15, phase.gasVolumeLPerKg(cellNo), ...
    partialPressure(phase.gasH2(cellNo)), partialPressure(phase.gasCO2(cellNo)), ...
    partialPressure(phase.gasCH4(cellNo)), partialPressure(phase.gasH2S(cellNo)), ...
    partialPressure(phase.gasN2(cellNo)), phase.pressureAtm(cellNo), ...
    phase.gasVolumeLPerKg(cellNo), phase.gasH2(cellNo)./phase.waterMass(cellNo), ...
    phase.gasCH4(cellNo)./phase.waterMass(cellNo), ...
    phase.gasCO2(cellNo)./phase.waterMass(cellNo), ...
    phase.gasH2S(cellNo)./phase.waterMass(cellNo), ...
    phase.gasN2(cellNo)./phase.waterMass(cellNo), ...
    mineral.calcite(cellNo), mineral.dolomite(cellNo), mineral.anhydrite(cellNo), ...
    mineral.quartz(cellNo), mineral.goethite(cellNo), mineral.pyrite(cellNo), ...
    mineral.brucite(cellNo), mineral.portlandite(cellNo), mineral.gypsum(cellNo), ...
    selectedOutputBlock());
end

function block = selectedOutputBlock()
% Deliberately no RATES or KINETICS: MRST is the reaction owner.
block = sprintf([ ...
    'SELECTED_OUTPUT 1\n' ...
    '-reset false\n' ...
    '-time true\n' ...
    '-step true\n' ...
    '-pH true\n' ...
    '-pe true\n' ...
    '-water true\n' ...
    '-molalities H2 N2 CarbonateO2 MethaneH4 H2Sulfide HCarbonateO3-\n' ...
    '-gases H2(g) redoxCarbonateO2(g) redoxCH4(g) redoxH2S(g) N2(g)\n' ...
    '-equilibrium_phases redoxCalcite redoxAnhydrite redoxGypsum redoxDolomite redoxGoethite redoxPyrite Brucite Portlandite Quartz\n' ...
    '-totals Carbonate(4) Sulfate(6) Ca Mg Fe_di Fe_tri K Na Cl Acetate Sulfide(-2) Si\n' ...
    '-saturation_indices redoxCalcite Brucite Portlandite redoxAnhydrite redoxDolomite Quartz redoxH2S(g)\n' ...
    'END\n']);
end

function raw = runIPhreeqcCOM(input, opt, cellNo)
try
    iph = actxserver(char(opt.comProgId));
catch ME
    error('H2Biochem:MrstMonodCOMActivation', ...
        'IPhreeqcCOM activation failed in cell %d for "%s":\n%s', ...
        cellNo, char(opt.comProgId), ME.message);
end
try
    loadStatus = iph.LoadDatabase(char(opt.databaseFile));
    if loadStatus ~= 0
        error('H2Biochem:MrstMonodCOMDatabase', ...
            'PHREEQC database load failed in cell %d:\n%s', cellNo, ...
            getPhreeqcError(iph));
    end
    status = iph.RunString(input);
    if status ~= 0
        error('H2Biochem:MrstMonodCOMRun', ...
            'PHREEQC equilibrium failed in cell %d:\n%s', cellNo, ...
            getPhreeqcError(iph));
    end
    raw = iph.GetSelectedOutputArray;
catch ME
    message = getPhreeqcError(iph);
    try
        clear iph
    catch
    end
    if startsWith(ME.identifier, 'H2Biochem:MrstMonodCOM')
        rethrow(ME);
    end
    error('H2Biochem:MrstMonodCOMRun', ...
        'IPhreeqcCOM failed in cell %d:\n%s\nPHREEQC message:\n%s', ...
        cellNo, ME.message, message);
end
clear iph
end

function message = getPhreeqcError(iph)
try
    message = char(iph.GetErrorString());
catch
    message = 'No PHREEQC error string was available from IPhreeqcCOM.';
end
if isempty(strtrim(message))
    message = 'IPhreeqcCOM did not provide an error message.';
end
end

function result = parseSelectedOutput(raw, cellNo)
context = sprintf('mrst-monod-com selected output in cell %d', cellNo);
[headers, values] = parseH2StorageIPhreeqcCOMSelectedOutput(raw, context);
read = @(aliases) getH2StorageIPhreeqcCOMSelectedOutputValue( ...
    headers, values, aliases, context);

result = emptyResult();
result.time = read({'time'});
result.step = read({'step'});
result.pH = read({'ph'});
result.pe = read({'pe'});
result.water = read({'massh2o', 'water'});
result.totalCarbon = read({'carbonate4molkgw', 'carbonate4'});
result.sulfate = read({'sulfate6molkgw', 'sulfate6'});
result.ca = read({'camolkgw', 'ca'});
result.mg = read({'mgmolkgw', 'mg'});
result.acetate = read({'acetatemolkgw', 'acetate'});
result.sulfide = read({'sulfide2molkgw', 'sulfide2'});
result.aqH2 = read({'mh2molkgw', 'mh2', 'h2'});
result.aqN2 = read({'mn2molkgw', 'mn2', 'n2'});
result.aqCO2 = read({'mcarbonateo2molkgw', 'mcarbonateo2', 'carbonateo2'});
result.aqCH4 = read({'mmethaneh4molkgw', 'mmethaneh4', 'methaneh4'});
result.aqH2S = read({'mh2sulfidemolkgw', 'mh2sulfide', 'h2sulfide'});
result.hco3 = read({'mhcarbonateo3molkgw', 'mhcarbonateo3', ...
    'mhco3molkgw', 'mhco3', 'hcarbonateo3', 'hco3'});
result.gasH2 = read({'gh2g', 'h2g'});
result.gasCO2 = read({'gredoxcarbonateo2g', 'redoxcarbonateo2g'});
result.gasCH4 = read({'gredoxch4g', 'redoxch4g'});
result.gasH2S = read({'gredoxh2sg', 'redoxh2sg'});
result.gasN2 = read({'gn2g', 'n2g'});
result.calcite = read({'redoxcalcite'});
result.anhydrite = read({'redoxanhydrite'});
result.gypsum = read({'redoxgypsum'});
result.dolomite = read({'redoxdolomite'});
result.goethite = read({'redoxgoethite'});
result.pyrite = read({'redoxpyrite'});
result.brucite = read({'brucite'});
result.portlandite = read({'portlandite'});
result.quartz = read({'quartz'});
end

function state = updateStateFromResults(model, state, opt, phase, result, waterMass)
nc = model.G.cells.num;
rhoWater = opt.waterDensity;
result = result(:);

state.phreeqcPH = column(result, 'pH');
state.phreeqcPE = column(result, 'pe');
state.phreeqcHCO3Molality = max(column(result, 'hco3'), 0);
state.phreeqcCO2Molality = max(column(result, 'aqCO2'), 0);
state.phreeqcTotalCarbon = max(column(result, 'totalCarbon'), 0);
state.tracerHCO3 = max(state.phreeqcTotalCarbon - state.phreeqcCO2Molality, 0).*rhoWater;
state.tracerSO4 = max(column(result, 'sulfate'), 0).*rhoWater;
state.tracerHS = max(column(result, 'sulfide'), 0).*rhoWater;
state.tracerCa = max(column(result, 'ca'), 0).*rhoWater;
state.tracerMg = max(column(result, 'mg'), 0).*rhoWater;

pka = state.phreeqcPH - log10(state.phreeqcHCO3Molality./ ...
    max(state.phreeqcCO2Molality, 1e-30));
previousPka = getStateVector(state, 'phreeqcCarbonatePka1', nc, ...
    model.carbonateBufferPka1);
pka(~isfinite(pka)) = previousPka(~isfinite(pka));
state.phreeqcCarbonatePka1 = pka;

state.phreeqcMineralCalcite = max(column(result, 'calcite'), 0).*waterMass;
state.phreeqcMineralDolomite = max(column(result, 'dolomite'), 0).*waterMass;
state.phreeqcMineralAnhydrite = max(column(result, 'anhydrite'), 0).*waterMass;
state.phreeqcMineralQuartz = max(column(result, 'quartz'), 0).*waterMass;
state.phreeqcMineralGoethite = max(column(result, 'goethite'), 0).*waterMass;
state.phreeqcMineralBrucite = max(column(result, 'brucite'), 0).*waterMass;
state.phreeqcMineralPortlandite = max(column(result, 'portlandite'), 0).*waterMass;
state.phreeqcMineralPyrite = max(column(result, 'pyrite'), 0).*waterMass;
state.phreeqcMineralGypsum = max(column(result, 'gypsum'), 0).*waterMass;

state.phreeqcMrstMonodComTime = column(result, 'time');
state.phreeqcMrstMonodComStep = column(result, 'step');
state.phreeqcMrstMonodComDissolvedH2Molality = max(column(result, 'aqH2'), 0);
state.phreeqcMrstMonodComDissolvedN2Molality = max(column(result, 'aqN2'), 0);
state.phreeqcMrstMonodComDissolvedCH4Molality = max(column(result, 'aqCH4'), 0);
state.phreeqcMrstMonodComDissolvedH2SMolality = max(column(result, 'aqH2S'), 0);
state.phreeqcMrstMonodComGasH2MolesPerKg = max(column(result, 'gasH2'), 0);
state.phreeqcMrstMonodComGasCO2MolesPerKg = max(column(result, 'gasCO2'), 0);
state.phreeqcMrstMonodComGasCH4MolesPerKg = max(column(result, 'gasCH4'), 0);
state.phreeqcMrstMonodComGasH2SMolesPerKg = max(column(result, 'gasH2S'), 0);
state.phreeqcMrstMonodComGasN2MolesPerKg = max(column(result, 'gasN2'), 0);

componentMoles = phase.componentMoles;
names = model.EOSModel.CompositionalMixture.names;
waterIndex = findComponent(names, {'H2O', 'Water'});
assert(~isempty(waterIndex), 'mrst-monod-com requires the EOS H2O component.');
waterMolarMass = model.EOSModel.CompositionalMixture.molarMass(waterIndex);
componentMoles(:, waterIndex) = max(componentMoles(:, waterIndex) + ...
    (column(result, 'water') - 1).*waterMass./waterMolarMass, 0);
componentMoles = replaceComponentMoles(componentMoles, names, {'H2', 'Hydrogen'}, ...
    (column(result, 'gasH2') + column(result, 'aqH2')).*waterMass, true);
componentMoles = replaceComponentMoles(componentMoles, names, {'CO2', 'CarbonDioxide'}, ...
    (column(result, 'gasCO2') + column(result, 'aqCO2')).*waterMass, true);
componentMoles = replaceComponentMoles(componentMoles, names, {'C1', 'CH4', 'Methane'}, ...
    (column(result, 'gasCH4') + column(result, 'aqCH4')).*waterMass, true);
componentMoles = replaceComponentMoles(componentMoles, names, {'H2S', 'HydrogenSulfide'}, ...
    (column(result, 'gasH2S') + column(result, 'aqH2S')).*waterMass, true);
componentMoles = replaceComponentMoles(componentMoles, names, {'N2', 'Nitrogen'}, ...
    (column(result, 'gasN2') + column(result, 'aqN2')).*waterMass, false);
componentMoles = replaceComponentMoles(componentMoles, names, ...
    {'CH3COOH', 'AceticAcid', 'Acetate'}, ...
    max(column(result, 'acetate'), 0).*waterMass, false);
assert(all(isfinite(componentMoles(:)) & componentMoles(:) >= 0), ...
    'mrst-monod-com returned invalid EOS component inventories.');

totalMoles = sum(componentMoles, 2);
assert(all(isfinite(totalMoles) & totalMoles > 0), ...
    'mrst-monod-com produced an invalid total EOS inventory.');
state.components = bsxfun(@rdivide, componentMoles, totalMoles);
state = clearStateFunctionCaches(model, state);
model = updateEOSSalinityForReflash(model, state);
state = model.computeFlash(state, inf);
state = refreshOutputStateFunctions(model, state);
state = updateDissolvedH2SLag(model, state);
end

function values = column(result, field)
values = reshape([result.(field)], [], 1);
end

function componentMoles = replaceComponentMoles(componentMoles, names, aliases, values, required)
index = findComponent(names, aliases);
if isempty(index)
    assert(~required, 'mrst-monod-com requires EOS component %s.', aliases{1});
    return;
end
componentMoles(:, index) = max(values, 0);
end

function componentMoles = getEOSComponentMoles(model, state)
nc = model.G.cells.num;
poreVolume = asCellVector(value(model.PVTPropertyFunctions.get( ...
    model, state, 'PoreVolume')), nc, 'pore volume');
[rhoL, rhoV, sL, sV] = getPhaseProperties(model, state);
totalMoles = poreVolume.*(sL.*rhoL + sV.*rhoV);
assert(all(isfinite(totalMoles) & totalMoles > 0), ...
    'mrst-monod-com requires positive EOS moles per cell.');
components = value(state.components);
assert(ismatrix(components) && size(components, 1) == nc && ...
    size(components, 2) == model.EOSModel.getNumberOfComponents(), ...
    'EOS components must be a cell-by-component matrix.');
componentMoles = bsxfun(@times, totalMoles, components);
end

function [rhoL, rhoV, sL, sV] = getPhaseProperties(model, state)
nc = model.G.cells.num;
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
propmodel = model.EOSModel.PropertyModel;
rhoL = propmodel.computeMolarDensity(model.EOSModel, value(state.pressure), ...
    value(state.x), value(state.Z_L), value(state.T), true);
rhoV = propmodel.computeMolarDensity(model.EOSModel, value(state.pressure), ...
    value(state.y), value(state.Z_V), value(state.T), false);
rhoL = asCellVector(value(rhoL), nc, 'liquid molar density');
rhoV = asCellVector(value(rhoV), nc, 'gas molar density');
sL = asCellVector(sL, nc, 'liquid saturation');
sV = asCellVector(sV, nc, 'gas saturation');
end

function values = getComponentMoles(componentMoles, names, aliases, required)
if nargin < 4
    required = true;
end
index = findComponent(names, aliases);
if isempty(index)
    assert(~required, 'mrst-monod-com requires EOS component %s.', aliases{1});
    values = zeros(size(componentMoles, 1), 1);
else
    values = componentMoles(:, index);
end
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
    model.EOSModel = model.EOSModel.enablesrb_coupling( ...
        value(state.tracerSO4), value(state.tracerHS), ...
        value(state.tracerHS) + value(state.h2sDissolvedLag), value(state.T));
end
end

function state = updateDissolvedH2SLag(model, state)
if ~model.sulfateReduction
    return;
end
names = model.EOSModel.CompositionalMixture.names;
index = find(strcmp(names, 'H2S'), 1);
if isempty(index)
    return;
end
x = value(state.x);
if iscell(x)
    xH2S = x{index};
else
    xH2S = x(:, index);
end
rhoL = model.EOSModel.PropertyModel.computeMolarDensity(model.EOSModel, ...
    value(state.pressure), value(state.x), value(state.Z_L), value(state.T), true);
state.h2sDissolvedLag = value(rhoL).*xH2S;
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

function values = getStateVector(state, field, nc, default)
if isfield(state, field)
    values = value(state.(field));
else
    values = default;
end
values = asCellVector(values, nc, field);
end

function values = optionVector(options, field, default, nc)
if isfield(options, field)
    values = options.(field);
else
    values = default;
end
values = asCellVector(values, nc, field);
end

function values = asCellVector(values, nc, name)
if isscalar(values)
    values = repmat(values, nc, 1);
else
    values = values(:);
    assert(numel(values) == nc, ...
        'mrst-monod-com expected %d %s values, got %d.', nc, name, numel(values));
end
assert(all(isfinite(values)), 'mrst-monod-com %s values must be finite.', name);
end

function index = findComponent(names, aliases)
index = [];
for i = 1:numel(aliases)
    index = find(strcmpi(names, aliases{i}), 1);
    if ~isempty(index)
        return;
    end
end
end

function validateSupportedComponents(names)
supported = {'H2O', 'Water', 'H2', 'Hydrogen', 'CO2', 'CarbonDioxide', ...
    'C1', 'CH4', 'Methane', 'H2S', 'HydrogenSulfide', 'N2', 'Nitrogen', ...
    'CH3COOH', 'AceticAcid', 'Acetate'};
unsupported = names(~cellfun(@(name) any(strcmpi(name, supported)), names));
assert(isempty(unsupported), ...
    ['mrst-monod-com cannot safely pass EOS components without a ', ...
     'PHREEQC_Modified.DAT mapping: %s'], strjoin(unsupported, ', '));
end

function isAbsolute = isAbsolutePath(path)
isAbsolute = ~isempty(regexp(path, '^[A-Za-z]:[\\/]|^\\\\', 'once')) || ...
    startsWith(path, filesep);
end

function result = emptyResult()
result = struct( ...
    'time', 0, 'step', 0, 'pH', 0, 'pe', 0, 'water', 0, ...
    'totalCarbon', 0, 'sulfate', 0, 'ca', 0, 'mg', 0, 'acetate', 0, ...
    'sulfide', 0, 'aqH2', 0, 'aqN2', 0, 'aqCO2', 0, 'aqCH4', 0, ...
    'aqH2S', 0, 'hco3', 0, 'gasH2', 0, 'gasCO2', 0, 'gasCH4', 0, ...
    'gasH2S', 0, 'gasN2', 0, 'calcite', 0, 'anhydrite', 0, ...
    'gypsum', 0, 'dolomite', 0, 'goethite', 0, 'pyrite', 0, ...
    'brucite', 0, 'portlandite', 0, 'quartz', 0);
end
