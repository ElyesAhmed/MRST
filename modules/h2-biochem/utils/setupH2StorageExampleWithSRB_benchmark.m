function [biochemFluid, model, schedule, state0] = setupH2StorageExampleWithSRB_benchmark(varargin)
% Set up the H2Storage1D benchmark with tracer-based sulfate reduction.
%
% This follows the H2Storage1D pure-H2 injection case, with either the
% moderate- or high-rate H2Storage1D kinetics.
% PHREEQC mineral reactions and pH evolution are excluded unless the
% optional post-timestep PHREEQC coupling is enabled.
%
% OPTIONAL PARAMETERS (property/value pairs):
%   rate                - 'medrate' (default) or 'highrate'. Selects the
%                         H2Storage1D metabolic kinetic-rate case; the
%                         well schedule is unchanged.
%   bacteriamodel       - Enable bacterial growth/decay (default: true)
%   bactDiffusion       - Enable microbial diffusion (default: false)
%   chemotaxisEffect    - Enable bacterial chemotaxis (default: false)
%   molecularDiffusion  - Enable molecular diffusion (default: false)
%   molecularDispersion - Enable mechanical dispersion (default: false)
%   bioClogging         - Enable bio-clogging (default: false)
%   carbonateBuffer     - Enable a finite fixed-pH HCO3-/CO2 buffer
%                         (default: false). Selects
%                         BiochemistryPhreeqcModel.
%   initialHCO3         - Initial bicarbonate molality when the buffer is
%                         enabled (default: Bentheimer value, 1.119e-3)
%   equilibrateInitialCO2 - Calibrate initial global CO2 so its flashed
%                         liquid mole fraction is in ideal fixed-pH
%                         equilibrium with initialHCO3 (default: false)
%   phreeqcDatabaseFile  - Absolute PHREEQC_Modified.DAT path used by
%                         either COM backend.
%   phreeqcBackend       - 'sequential-compositional-phreeqc' (default,
%                         Windows IPhreeqcCOM kinetics), or
%                         'sequential-h2biochem-phreeqc' (Windows IPhreeqcCOM equilibrium).
%   phreeqcComProgId     - Registered Windows IPhreeqcCOM ProgID for
%                         COM backends (default: IPhreeqcCOM.Object)
%   phreeqc*WtFraction   - Dolomite, calcite, brucite, and quartz rock
%                         weight fractions used to seed PHREEQC minerals
%   phreeqcTimestepCoupling - After every converged timestep, equilibrate
%                         each cell with the selected PHREEQC backend
%                         (default: false).
%                         sequential-h2biochem-phreeqc must be run with
%                         simulateSequentialH2BiochemPhreeqc, not direct
%                         simulateScheduleAD.
%   phreeqcPicardMaxIterations, phreeqcPicardRelaxation,
%   phreeqcPicardAbsoluteTolerance, phreeqcPicardRelativeTolerance,
%   phreeqcPicardPHTolerance, phreeqcPicardPkaTolerance - sequential-h2biochem-phreeqc
%                         outer Picard controls.
%   phreeqcElementBalanceAbsoluteTolerance,
%   phreeqcElementBalanceRelativeTolerance - Scalar or six-element
%                         H/C/S/Ca/Mg/Fe conservation tolerances.
%   sequentialCompositionalPhreeqcInitialC4 - Initial nonvolatile C(4)
%                         molality (default: 1.370e-3).
%   paperBiomassKinetics - Use first-order biomass decay and the paper's
%                         Nmax/N0 = 1e4 population range (default: false).
%                         Selects BiochemistryPhreeqcModel because this
%                         benchmark-specific kinetics option is not part of
%                         the original BiochemistryModel.
%   nbact0              - Initial normalized bacterial state, scalar or
%                         [MET ACE SRB] vector (default: 60.6)
%   injectionCO2        - CO2 mole fraction in injected gas (default: 0,
%                         matching the paper; nonzero values are extensions)

require ad-props compositional deckformat h2-biochem

opt = struct(...
    'rate', 'highrate', ...
    'bacteriamodel', true, ...
    'bactDiffusion', false, ...
    'chemotaxisEffect', false, ...
    'molecularDiffusion', false, ...
    'molecularDispersion', false, ...
    'bioClogging', false, ...
    'carbonateBuffer', false, ...
    'initialHCO3', 1.1119e-3, ...
    'equilibrateInitialCO2', false, ...
    'phreeqcTimestepCoupling', false, ...
    'phreeqcBackend', 'sequential-compositional-phreeqc', ...
    'phreeqcDatabaseFile', '', ...
    'phreeqcComProgId', 'IPhreeqcCOM.Object', ...
    'phreeqcPicardMaxIterations', 1, ...
    'phreeqcPicardRelaxation', 0.9, ...
    'phreeqcPicardAbsoluteTolerance', 1e-8, ...
    'phreeqcPicardRelativeTolerance', 1e-2, ...
    'phreeqcPicardPHTolerance', 2e-2, ...
    'phreeqcPicardPkaTolerance', 2e-2, ...
    'sequentialCompositionalPhreeqcInitialC4', 1.370e-3, ...
    'phreeqcPicardReactionAbsoluteTolerance', 1e-12, ...
    'phreeqcPicardReactionRelativeTolerance', 1e-3, ...
    'phreeqcElementBalanceAbsoluteTolerance', 1e-7, ...
    'phreeqcElementBalanceRelativeTolerance', 1e-8, ...
    'phreeqcDolomiteWtFraction', 0.02, ...
    'phreeqcCalciteWtFraction', 0, ...
    'phreeqcBruciteWtFraction', 0, ...
    'phreeqcQuartzWtFraction', 0.98, ...
    'paperBiomassKinetics', false, ...
    'nbact0', 60.6, ...
    'injectionCO2', 0.0);
opt = merge_options(opt, varargin{:});

assert(opt.injectionCO2 >= 0 && opt.injectionCO2 <= 1, ...
    'injectionCO2 must be a mole fraction between zero and one.');
assert(ischar(opt.phreeqcBackend) || ...
    (isstring(opt.phreeqcBackend) && isscalar(opt.phreeqcBackend)), ...
    'phreeqcBackend must be a character vector or scalar string.');
phreeqcBackend = lower(strtrim(char(opt.phreeqcBackend)));
assert(ismember(phreeqcBackend, ...
    {'sequential-compositional-phreeqc', 'sequential-h2biochem-phreeqc'}), ...
    ['phreeqcBackend must be ''sequential-compositional-phreeqc'' or ', ...
     '''sequential-h2biochem-phreeqc''.']);
if opt.phreeqcTimestepCoupling
    assert(opt.carbonateBuffer && opt.bacteriamodel, ...
        ['phreeqcTimestepCoupling requires carbonateBuffer=true and ', ...
         'bacteriamodel=true so HCO3 and SO4 tracers are active.']);
    assert(ispc, ['COM PHREEQC backends require Windows and ', ...
        'a registered IPhreeqcCOM server.']);
    assert(~isempty(opt.phreeqcDatabaseFile) && ...
        isAbsolutePhreeqcPath(opt.phreeqcDatabaseFile), ...
        ['COM PHREEQC backends require phreeqcDatabaseFile to be an explicit ', ...
         'absolute path to PHREEQC_Modified.DAT.']);
    assert(isfile(opt.phreeqcDatabaseFile), ...
        'COM PHREEQC database not found: %s', opt.phreeqcDatabaseFile);
    assert(contains(lower(char(opt.phreeqcDatabaseFile)), 'phreeqc_modified.dat'), ...
        ['COM PHREEQC backends require PHREEQC_Modified.DAT, not a ', ...
         'standard PHREEQC database: %s'], opt.phreeqcDatabaseFile);
    balanceTolerances = {'phreeqcElementBalanceAbsoluteTolerance', ...
        'phreeqcElementBalanceRelativeTolerance'};
    for i = 1:numel(balanceTolerances)
        name = balanceTolerances{i};
        values = opt.(name);
        validateattributes(values, {'numeric'}, ...
            {'vector', 'real', 'finite', 'nonnegative'}, mfilename, name);
        assert(isscalar(values) || numel(values) == 6, ...
            '%s must be scalar or contain H/C/S/Ca/Mg/Fe values.', name);
    end
end
if strcmp(phreeqcBackend, 'sequential-h2biochem-phreeqc')
    validateattributes(opt.phreeqcPicardMaxIterations, {'numeric'}, ...
        {'scalar', 'integer', 'finite', '>=', 1}, mfilename, ...
        'phreeqcPicardMaxIterations');
    validateattributes(opt.phreeqcPicardRelaxation, {'numeric'}, ...
        {'scalar', 'real', 'finite', '>', 0, '<=', 1}, mfilename, ...
        'phreeqcPicardRelaxation');
    picardTolerances = {'phreeqcPicardAbsoluteTolerance', ...
        'phreeqcPicardRelativeTolerance', 'phreeqcPicardPHTolerance', ...
        'phreeqcPicardReactionAbsoluteTolerance', ...
        'phreeqcPicardReactionRelativeTolerance'};
    for i = 1:numel(picardTolerances)
        validateattributes(opt.(picardTolerances{i}), {'numeric'}, ...
            {'scalar', 'real', 'finite', 'positive'}, mfilename, picardTolerances{i});
    end
end

rateCase = lower(char(opt.rate));
switch rateCase
    case {'medrate', 'moderate'}
        metabolicGrowthRate = [1.109, 0.872, 1.048]./86400;
        h2HalfSaturation = [10e-6, 2.5e-6, 2.9e-6]./55.5;
        rateCase = 'medrate';
    case {'highrate', 'high'}
        metabolicGrowthRate = [4.1, 1.9, 5.5]./86400;
        h2HalfSaturation = [9e-6, 2.5e-6, 2.9e-6]./55.5;
        rateCase = 'highrate';
    otherwise
        error('rate must be ''medrate'' or ''highrate''.');
end

%% H2Storage1D grid, rock, and transport properties
G = computeGeometry(cartGrid([50, 1, 1], [50, 1, 1]));
rock = makeRock(G, 100*milli*darcy, 0.20);

Sw_table = [0, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, ...
    0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 0.999];
krW_table = [0, 0, 0.002, 0.01, 0.02, 0.033, 0.049, 0.066, 0.09, 0.119, 0.15, ...
    0.186, 0.227, 0.277, 0.33, 0.39, 0.462, 0.54, 0.62, 0.71, 0.8, 0.9, 1];
Sg_table = [0.05, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, ...
    0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 1];
krG_table = [0, 0.013, 0.026, 0.04, 0.058, 0.078, 0.1, 0.126, 0.156, 0.187, ...
    0.222, 0.26, 0.3, 0.348, 0.4, 0.45, 0.505, 0.562, 0.62, 0.68, 0.74, 1];

fluid = initSimpleADIFluid('phases', 'OG', ...
    'mu', [1.3059*centi*poise, 0.01763*centi*poise], ...
    'rho', [999.7, 1.2243].*kilogram/meter^3, ...
    'pRef', 150*barsa, ...
    'c', [5e-5/barsa, 0.0067/barsa], ...
    'n', [2, 2], ...
    'smin', [0.2, 0.05]);
fluid.krO = @(sw) interpTable(Sw_table, krW_table, sw);
fluid.krG = @(sg) interpTable(Sg_table, krG_table, sg);
fluid.pcOG = @(sg) 0;

%% Components, reactions, and fixed brine chemistry
volatileNames = {'Water', 'Hydrogen', 'CarbonDioxide', 'Methane', ...
    'HydrogenSulfide', 'AceticAcid'};
volatileSymbols = {'H2O', 'H2', 'CO2', 'C1', 'H2S', 'CH3COOH'};
compFluid = TableCompositionalMixture(volatileNames, volatileSymbols);

reactNames = {'MethanogenicArchae', 'AcetogenicBacteria', ...
    'SulfateReducingBacteria'};
biomassNames = {'bactM', 'bactA', 'bactS'};
biochemFluid = TableBioChemMixture(reactNames, biomassNames);

% H2Storage1D kinetics. Yield scales and bacterial-state limits retain the
% h2-biochem formulation's calibration in both rate cases.
idxM = strcmp(biochemFluid.metabolicReaction, 'MethanogenicArchae');
biochemFluid.Psigrowthmax(idxM) = metabolicGrowthRate(1);
biochemFluid.alphaH2(idxM) = h2HalfSaturation(1);
biochemFluid.alphasub(idxM) = 230e-6/55.5;
biochemFluid.bbact(idxM) = 0.01*biochemFluid.Psigrowthmax(idxM);
biochemFluid.Y_H2(idxM) = 0.03*3.333e12;
biochemFluid.nbactMax(idxM) = 1e10;

idxA = strcmp(biochemFluid.metabolicReaction, 'AcetogenicBacteria');
biochemFluid.Psigrowthmax(idxA) = metabolicGrowthRate(2);
biochemFluid.alphaH2(idxA) = h2HalfSaturation(2);
biochemFluid.alphasub(idxA) = 115.5e-6/55.5;
biochemFluid.bbact(idxA) = 0.01*biochemFluid.Psigrowthmax(idxA);
biochemFluid.Y_H2(idxA) = 0.07*3.333e12;
biochemFluid.nbactMax(idxA) = 1e10;

idxS = strcmp(biochemFluid.metabolicReaction, 'SulfateReducingBacteria');
biochemFluid.gamrH2(idxS) = -4;
biochemFluid.gamrsub(idxS) = -1;
biochemFluid.gampH2O(idxS) = 4;
biochemFluid.gamp2(idxS) = 1;
biochemFluid.Psigrowthmax(idxS) = metabolicGrowthRate(3);
biochemFluid.alphaH2(idxS) = h2HalfSaturation(3);
biochemFluid.alphasub(idxS) = 2751.5e-6/55.5;
biochemFluid.bbact(idxS) = 0.01*biochemFluid.Psigrowthmax(idxS);
biochemFluid.Y_H2(idxS) = 0.08*3.333e12;
biochemFluid.nbactMax(idxS) = 1e9;

% H2Storage1D initial pH, Na concentration used as NaCl-equivalent
% salinity, and S6 sulfate concentration. These remain fixed except for
% sulfate/HS tracer transport and reaction; no PHREEQC chemistry is used.
pH0 = 6.24;
initialNaCl = 2.865;
initialSO4 = 4.664e-3;
initialHCO3 = opt.initialHCO3;
if opt.phreeqcTimestepCoupling
    % Coupled cases initialize Solution.C4 to 1.370e-3 mol/kgw.
    validateattributes(opt.sequentialCompositionalPhreeqcInitialC4, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, mfilename, 'sequentialCompositionalPhreeqcInitialC4');
    initialHCO3 = opt.sequentialCompositionalPhreeqcInitialC4;
end
carbonateBufferPka1 = 6.35;
eos = SoreideWhitsonEos(G, compFluid, ...
    'msalt', initialNaCl, ...
    'pH', pH0, ...
    'initial_NaCl', initialNaCl, ...
    'initial_SO4', initialSO4, ...
    'rho_water', 1000);
phreeqcCouplingOptions = struct( ...
    'databaseFile', opt.phreeqcDatabaseFile, ...
    'comProgId', opt.phreeqcComProgId, ...
    'waterDensity', eos.rho_water, ...
    'phreeqcElementBalanceAbsoluteTolerance', ...
        opt.phreeqcElementBalanceAbsoluteTolerance, ...
    'phreeqcElementBalanceRelativeTolerance', ...
        opt.phreeqcElementBalanceRelativeTolerance, ...
    'Na', 2.865, ...
    'Ca', 0.2857, ...
    'Mg', 0.1144, ...
    'Cl', 3.655, ...
    'Si', 9.723e-5, ...
    'dolomitePhase', 'Dolomite', ...
    'calcitePhase', 'Calcite', ...
    'brucitePhase', 'Brucite', ...
    'quartzPhase', 'Quartz', ...
    'dolomiteWtFraction', opt.phreeqcDolomiteWtFraction, ...
    'calciteWtFraction', opt.phreeqcCalciteWtFraction, ...
    'bruciteWtFraction', opt.phreeqcBruciteWtFraction, ...
    'quartzWtFraction', opt.phreeqcQuartzWtFraction, ...
    'initialDolomiteMoles', [], ...
    'initialCalciteMoles', [], ...
    'initialBruciteMoles', [], ...
    'initialQuartzMoles', [], ...
    'initialDolomiteMolality', [], ...
    'initialCalciteMolality', [], ...
    'initialBruciteMolality', [], ...
    'initialQuartzMolality', []);
phreeqcCouplingOptions.sequentialCompositionalPhreeqcMuMET = metabolicGrowthRate(1)*day;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcMuACE = metabolicGrowthRate(2)*day;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcMuSRB = metabolicGrowthRate(3)*day;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcBMET = 0.01*phreeqcCouplingOptions.sequentialCompositionalPhreeqcMuMET;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcBACE = 0.01*phreeqcCouplingOptions.sequentialCompositionalPhreeqcMuACE;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcBSRB = 0.01*phreeqcCouplingOptions.sequentialCompositionalPhreeqcMuSRB;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcYMET = 0.03*4;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcYACE = 0.07*4;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcYSRB = 0.08*4;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcKDMET = h2HalfSaturation(1)*55.5;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcKDACE = h2HalfSaturation(2)*55.5;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcKDSRB = h2HalfSaturation(3)*55.5;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcKAMET = 230e-6;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcKAACE = 115.5e-6;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcKASRB = 2751.5e-6;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcN0 = 1e9;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcNmax = 1e13;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcCellMass = 1e-14;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcBiomassMW = 24.6;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcSteps = 5;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcInitialAnhydriteMoles = 0;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcInitialGoethiteMoles = 0;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcInitialPortlanditeMoles = 0;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcInitialPyriteMoles = 0;
phreeqcCouplingOptions.sequentialCompositionalPhreeqcInitialGypsumMoles = 0;
phreeqcCouplingOptions.sequentialH2BiochemPhreeqcMaxIterations = opt.phreeqcPicardMaxIterations;
phreeqcCouplingOptions.sequentialH2BiochemPhreeqcRelaxation = opt.phreeqcPicardRelaxation;
phreeqcCouplingOptions.sequentialH2BiochemPhreeqcAbsoluteTolerance = ...
    opt.phreeqcPicardAbsoluteTolerance;
phreeqcCouplingOptions.sequentialH2BiochemPhreeqcRelativeTolerance = ...
    opt.phreeqcPicardRelativeTolerance;
phreeqcCouplingOptions.sequentialH2BiochemPhreeqcPHTolerance = opt.phreeqcPicardPHTolerance;
phreeqcCouplingOptions.sequentialH2BiochemPhreeqcPkaTolerance = opt.phreeqcPicardPkaTolerance;
phreeqcCouplingOptions.sequentialH2BiochemPhreeqcReactionAbsoluteTolerance = ...
    opt.phreeqcPicardReactionAbsoluteTolerance;
phreeqcCouplingOptions.sequentialH2BiochemPhreeqcReactionRelativeTolerance = ...
    opt.phreeqcPicardReactionRelativeTolerance;
%% Model assembly
backend = DiagonalAutoDiffBackend('modifyOperators', true);
modelOptions = { ...
    'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', opt.bacteriamodel, ...
    'bactDiffusion', opt.bactDiffusion, ...
    'chemotaxisEffect', opt.chemotaxisEffect, ...
    'molecularDiffusion', opt.molecularDiffusion, ...
    'molecularDispersion', opt.molecularDispersion, ...
    'enableSulfateSource', false, ...
    'liquidPhase', 'O', 'vaporPhase', 'G'};
usePhreeqcModel = opt.phreeqcTimestepCoupling || opt.carbonateBuffer || ...
    opt.paperBiomassKinetics;
if usePhreeqcModel
    model = BiochemistryPhreeqcModel(G, rock, fluid, compFluid, ...
        biochemFluid, true, backend, modelOptions{:}, ...
        'carbonateBuffer', opt.carbonateBuffer, ...
        'carbonateBufferPH', pH0, ...
        'carbonateBufferPka1', carbonateBufferPka1, ...
        'phreeqcTimestepCoupling', opt.phreeqcTimestepCoupling, ...
        'phreeqcBackend', phreeqcBackend, ...
        'phreeqcDatabaseFile', opt.phreeqcDatabaseFile, ...
        'phreeqcComProgId', opt.phreeqcComProgId, ...
        'phreeqcCouplingOptions', phreeqcCouplingOptions, ...
        'bacterialDecayOrder', 1 + ~opt.paperBiomassKinetics);
else
    model = BiochemistryModel(G, rock, fluid, compFluid, biochemFluid, ...
        true, backend, modelOptions{:});
end
model.EOSModel = eos;
model.OutputStateFunctions{end + 1} = 'ComponentPhaseDensity';
if opt.bacteriamodel
    model.OutputStateFunctions{end + 1} = 'PsiGrowthRate';
    model.OutputStateFunctions{end + 1} = 'BacterialMass';
    if usePhreeqcModel
        model.OutputStateFunctions{end + 1} = 'CarbonLimitedGrowthRate';
    end
end
nBioReactions = biochemFluid.nbioreact;
validateattributes(opt.nbact0, {'numeric'}, ...
    {'real', 'finite', 'positive', 'vector'}, mfilename, 'nbact0');
if isscalar(opt.nbact0)
    nbact0 = repmat(opt.nbact0, 1, nBioReactions);
else
    assert(numel(opt.nbact0) == nBioReactions, ...
        'nbact0 must be scalar or contain one value per biochemical reaction.');
    nbact0 = reshape(opt.nbact0, 1, []);
end
if opt.bioClogging && opt.bacteriamodel
    model = setupBioCloggingModel(model, nbact0, [180, 180, 180], [0.5, 0.5, 0.5], true);
else
    model = setupBioCloggingModel(model, nbact0, [180, 180, 180], [0, 0, 0], false);
end

% H2Storage1D starts with zero anhydrite. Disable the h2-biochem
% anhydrite-dissolution extension so the sulfate tracer is not replenished.
model = model.validateModel();
if opt.bacteriamodel
    model = setBenchmarkSulfateSource(model, false);
end
%% H2Storage1D initial fluid state
p0 = 150*barsa;
T0 = 273.15 + 60;
zCO2 = 0.0025;
if opt.phreeqcTimestepCoupling
    % H2Storage1D starts with no EOS CO2; C(4) is represented solely by
    % Solution.C4 and must not be duplicated in the volatile inventory.
    zCO2 = 0;
elseif opt.equilibrateInitialCO2
    assert(opt.carbonateBuffer && opt.initialHCO3(1) > 0, ...
        ['equilibrateInitialCO2 requires carbonateBuffer=true and ', ...
         'a positive initialHCO3.']);
    [zCO2, xCO2Target] = calibrateInitialCO2(model, p0, T0, nbact0, eos, ...
        opt.initialHCO3./10.^(pH0 - carbonateBufferPka1), opt.bacteriamodel);
end
z0 = zeros(model.G.cells.num, compFluid.getNumberOfComponents());
z0(:,strcmp(compFluid.names, 'H2O')) = 0.90;
z0(:,strcmp(compFluid.names, 'C1')) = 0.10 - zCO2;
z0(:,strcmp(compFluid.names, 'CO2')) = zCO2;
if opt.bacteriamodel
    state0 = initCompositionalStateBacteria(model, p0, T0, [], z0, nbact0, eos);
else
    state0 = initCompositionalState(model, p0, T0, [], z0, eos);
end

if opt.bacteriamodel
    ncell = G.cells.num;
    state0.tracerSO4 = repmat(initialSO4*eos.rho_water, ncell, 1);
    state0.tracerHS = zeros(ncell, 1);
    state0.h2sDissolvedLag = zeros(ncell, 1);
    if usePhreeqcModel && opt.carbonateBuffer
        initialNonvolatileDIC = initialHCO3;
        state0.tracerHCO3 = repmat(max(initialNonvolatileDIC, 0)*eos.rho_water, ncell, 1);
    end
end
% Add extra initial aqueous carbon
% extraC4 = e-1; % mol/kg water
% state0.tracerHCO3 = [0.5550
%     0.5395
%     0.5286
%     0.5225
%     0.5193
%     0.5176
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701
%     1.3701]; 
% state0.tracerHCO3= state0.tracerHCO3 + ...
%     extraC4*model.EOSModel.rho_water;
if opt.phreeqcTimestepCoupling
    initialCa = phreeqcCouplingOptions.Ca;
    initialMg = phreeqcCouplingOptions.Mg;
    state0.tracerCa = repmat(initialCa*eos.rho_water, ncell, 1);
    state0.tracerMg = repmat(initialMg*eos.rho_water, ncell, 1);
    state0.phreeqcHCO3Molality = repmat(initialHCO3, ncell, 1);
    state0.phreeqcPE = 4*ones(ncell, 1);
    state0 = initializeH2StoragePhreeqcCouplingState(model, state0);
end
%% H2Storage1D benchmark: 50 d injection, 150 d storage, 50 d production.
dtInj = repmat(2*day, 25, 1);
dtStor = repmat(2*day, 75, 1);
dtProd = repmat(2*day, 25, 1);

schedule.step.val = [dtInj; dtStor; dtProd];
schedule.step.control = [ ...
    ones(numel(dtInj), 1); ...
    2*ones(numel(dtStor), 1); ...
    3*ones(numel(dtProd), 1)];

[~, ~, ~, ~, Z_V] = standaloneFlash(p0, T0, [0, 1, 0, 0, 0, 0], eos);
Bg = 101325/298.15*Z_V*T0/p0;
pv = sum(G.cells.volumes.*rock.poro);
rate = 0.005*pv*meter^3/day/Bg;

injComponents = [0, 1 - opt.injectionCO2, opt.injectionCO2, 0, 0, 0];
W1 = verticalWell([], G, rock, 1, 1, 1, ...
    'Type', 'rate', 'Val', rate, 'Name', 'Injector', ...
    'comp_i', [0, 1], 'sign', 1, 'radius', 0.1);
W1(1).components = injComponents;

W2 = verticalWell([], G, rock, 1, 1, 1, ...
    'Type', 'rate', 'Val', 0, 'Name', 'Shut-in', ...
    'comp_i', [0, 1], 'sign', 1, 'radius', 0.1);
W2(1).components = [0, 1, 0, 0, 0, 0];

W3 = verticalWell([], G, rock, 1, 1, 1, ...
    'Type', 'grat', 'Val', -rate, 'Name', 'Producer', ...
    'comp_i', [0.5, 0.5], 'sign', -1, 'radius', 0.1);
W3(1).components = [0, 1, 0, 0, 0, 0];
W3.lims.bhp = p0;

tmp = cell(4, 1);
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);
schedule.control(1).W = W1;
schedule.control(2).W = W2;
schedule.control(3).W = W3;

fprintf('H2Storage1D SRB tracer benchmark: pH %.2f, SO4 %.4g mol/kgw, NaCl-equivalent %.3f mol/kgw\n', ...
    pH0, initialSO4, initialNaCl);
fprintf('Kinetic rate case: %s\n', rateCase);
if usePhreeqcModel && opt.carbonateBuffer
    fprintf('Carbonate buffer: %.4g mol/kgw HCO3- at pH %.2f\n', ...
        initialHCO3, pH0);
end
if ~opt.phreeqcTimestepCoupling && opt.equilibrateInitialCO2
    fprintf(['Initial CO2 calibrated to x_CO2(liquid) %.4g: ', ...
        'z_CO2 %.4g\n'], xCO2Target, zCO2);
end
if opt.phreeqcTimestepCoupling
    fprintf(['PHREEQC timestep coupling (%s): post-convergence, per-cell ', ...
        'chemistry split\n'], phreeqcBackend);
    if strcmp(phreeqcBackend, 'sequential-compositional-phreeqc')
        fprintf(['compositional PHREEQC COM kinetics use separate PHREEQC biomass state; ', ...
            'implicit h2-biochem reaction sources are disabled.\n']);
    elseif strcmp(phreeqcBackend, 'sequential-h2biochem-phreeqc')
        fprintf(['MRST nbact Monod sources remain the sole reaction owner. ', ...
            'Run simulateSequentialH2BiochemPhreeqc for same-timestep ', ...
            'equilibrium feedback (max %d Picard iterations).\n'], ...
            opt.phreeqcPicardMaxIterations);
    end
end
fprintf('Injection: %.1f%% H2 + %.1f%% CO2\n', ...
    100*(1 - opt.injectionCO2), 100*opt.injectionCO2);
end

function [zCO2, xCO2Target] = calibrateInitialCO2(model, p0, T0, nbact0, eos, mCO2, useBacteria)
% mCO2 : scalar or vector of target aqueous CO2 molality [mol/kgw].
% Returns zCO2 (gas mole fraction) and xCO2Target (liquid mole fraction) of same size.

iCO2 = find(strcmp(model.EOSModel.CompositionalMixture.names, 'CO2'), 1);
xCO2Target = mCO2 ./ (55.508 + mCO2);   % vector if mCO2 is vector

lower = 0;
upper = 0.10 - 1e-10;

% Pre-allocate output
zCO2 = zeros(size(xCO2Target));

% Residual function that takes candidate and target (scalar)
    function res = residual(candidate, target)
        z = zeros(1, numel(model.EOSModel.CompositionalMixture.names));
        z(strcmp(model.EOSModel.CompositionalMixture.names, 'H2O')) = 0.90;
        z(strcmp(model.EOSModel.CompositionalMixture.names, 'C1')) = 0.10 - candidate;
        z(iCO2) = candidate;
        if useBacteria
            state = initCompositionalStateBacteria(model, p0, T0, [], z, nbact0, eos);
        else
            state = initCompositionalState(model, p0, T0, [], z, eos);
        end
        res = state.x(1, iCO2) - target;   % scalar
    end

% Bracket check for each target (element-wise)
fLower = arrayfun(@(t) residual(lower, t), xCO2Target);
fUpper = arrayfun(@(t) residual(upper, t), xCO2Target);
assert(all(fLower .* fUpper <= 0), ...
    'Bracket check failed for at least one target. No valid zCO2 found.');

% Solve for each target – using a loop for clarity and robustness
for i = 1:numel(xCO2Target)
    target = xCO2Target(i);
    % Bind the target into a scalar residual for fzero
    zCO2(i) = fzero(@(c) residual(c, target), [lower, upper]);
end
end


function model = setBenchmarkSulfateSource(model, enabled)
model.enableSulfateSource = enabled;
tracerRate = model.FlowDiscretization.getStateFunction('SRBTracerConvRate');
tracerRate.enable_sulfate_source = enabled;
model.FlowDiscretization = model.FlowDiscretization.setStateFunction( ...
    'SRBTracerConvRate', tracerRate);

facilityRate = model.FacilityModel.FacilityFlowDiscretization.getStateFunction( ...
    'SRBTracerConvRate');
facilityRate.enable_sulfate_source = enabled;
model.FacilityModel.FacilityFlowDiscretization = ...
    model.FacilityModel.FacilityFlowDiscretization.setStateFunction( ...
        'SRBTracerConvRate', facilityRate);
end

function isAbsolute = isAbsolutePhreeqcPath(path)
path = char(path);
isAbsolute = ~isempty(regexp(path, '^[A-Za-z]:[\\/]|^\\\\', 'once')) || ...
    startsWith(path, filesep);
end
