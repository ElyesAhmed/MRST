function [biochemFluid, model, schedule, state0] = setupH2StorageExampleWithSRB_benchmark(varargin)
% Set up the H2Storage1D benchmark with tracer-based sulfate reduction.
%
% This retains the h2-biochem reactions and the 95/5 H2/CO2 injection
% extension, with either the moderate- or high-rate H2Storage1D kinetics.
% PHREEQC mineral reactions and pH evolution are intentionally excluded.
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
%   nbact0              - Initial normalized bacterial state (default: 60.6)
%   injectionCO2        - CO2 mole fraction in injected gas (default: 0.05)

require ad-props compositional deckformat h2-biochem

opt = struct(...
    'rate', 'highrate', ...
    'bacteriamodel', true, ...
    'bactDiffusion', false, ...
    'chemotaxisEffect', false, ...
    'molecularDiffusion', false, ...
    'molecularDispersion', false, ...
    'bioClogging', false, ...
    'nbact0', 60.6, ...
    'injectionCO2', 0.025);
opt = merge_options(opt, varargin{:});

assert(opt.injectionCO2 >= 0 && opt.injectionCO2 <= 1, ...
    'injectionCO2 must be a mole fraction between zero and one.');

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
eos = SoreideWhitsonEos(G, compFluid, ...
    'msalt', initialNaCl, ...
    'pH', pH0, ...
    'initial_NaCl', initialNaCl, ...
    'initial_SO4', initialSO4, ...
    'rho_water', 1000);

%% Model assembly
backend = DiagonalAutoDiffBackend('modifyOperators', true);
model = BiochemistryModel(G, rock, fluid, compFluid, biochemFluid, true, backend, ...
    'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', opt.bacteriamodel, ...
    'bactDiffusion', opt.bactDiffusion, ...
    'chemotaxisEffect', opt.chemotaxisEffect, ...
    'molecularDiffusion', opt.molecularDiffusion, ...
    'molecularDispersion', opt.molecularDispersion, ...
    'enableSulfateSource', false, ...
    'liquidPhase', 'O', 'vaporPhase', 'G');
model.EOSModel = eos;
model.OutputStateFunctions{end + 1} = 'ComponentPhaseDensity';
if opt.bacteriamodel
    model.OutputStateFunctions{end + 1} = 'PsiGrowthRate';
    model.OutputStateFunctions{end + 1} = 'BacterialMass';
end
nBioReactions = biochemFluid.nbioreact;
nbact0 = opt.nbact0*ones(1, nBioReactions);
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
z0 = zeros(1, compFluid.getNumberOfComponents());
z0(strcmp(compFluid.names, 'H2O')) = 0.90;
z0(strcmp(compFluid.names, 'C1')) = 0.10;

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
end
%% H2Storage1D schedule with the retained H2/CO2 injection extension
dt = repmat(2*day, 125, 1);
schedule = simpleSchedule(dt);
schedule.step.control(26:100) = 2;
schedule.step.control(101:end) = 3;

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
fprintf('Injection: %.1f%% H2 + %.1f%% CO2\n', ...
    100*(1 - opt.injectionCO2), 100*opt.injectionCO2);
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
