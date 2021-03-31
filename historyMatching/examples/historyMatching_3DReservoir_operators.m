%% History matching of 3D reservoir simulation
% This is a minimal example (and a sanity test) for creating an ensemble of
% 3D reservoirs with two phases (oil and water), when the uncertain
% parameters are transmissibility and porevolume defined through the
% OperatorSamples class
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble 

mrstVerbose off

%% Set up and simulate the true solution
% We will here use an identical twin experiment, where we use the same
% problem for both generating the truth and as a base for our ensemble.

trueProblemName = 'ensemble_base_problem_3d_reservoir';
baseProblemOptions = {};


directoryTruth = fullfile(mrstOutputDirectory(), ...
                          'historyMatching', 'truth', 'operators', ...
                          trueProblemName);
                      
topDirectory = fullfile(mrstOutputDirectory(), ...
                        'historyMatching', 'tutorial', 'operators', trueProblemName);
                      
trueExample = MRSTExample(trueProblemName);
trueProblem = trueExample.getPackedSimulationProblem('Directory', directoryTruth);

plotExample = false;
rerunTrueProblemFromScratch = false;


if rerunTrueProblemFromScratch
    clearPackedSimulatorOutput(trueProblem);
end
simulatePackedProblem(trueProblem);
if plotExample
    [wellSols, states, reports] = getPackedSimulatorOutput(trueProblem);
    trueExample.plot(states);
end

%% Generate observations
% Define a QoI object for storing the relevant observations we will use for
% history matching

trueQoI = WellQoIHM(...
    'wellNames', {'P1', 'P2'}, ...
    'fldname', {'qOs', 'qWs'}, ...
    'cumulative', false);

trueQoI = trueQoI.validateQoI(trueProblem);
trueObservations = trueQoI.getQoI(trueProblem);

% Define observation uncertainty 
obsStdDev = 0.0004*0.1;

% Create a separate ResultHandler for the observations 
observationResultHandler = trueQoI.ResultHandler;
observationResultHandler.dataPrefix = 'observedQoI';

% Add some observation noise and store output
if numel(observationResultHandler.getValidIds) < 1
    for w = 1:numel(trueQoI.wellNames)
        for f = 1:numel(trueQoI.fldname)
            perturbedObservations{w}{f} = trueObservations{w}{f} + randn(size(trueObservations{w}{f}))*obsStdDev;
        end
    end
    observationResultHandler{1} = {perturbedObservations};
end



%% Select and populate samples for stochastic configurations class

ensembleSize = 70;

% We model the transmisibility as a log-Gaussian distribution with the same
% mean and standard deviation as in trueExample

logMean = mean(log(trueExample.model.operators.T));
logStd = std(log(trueExample.model.operators.T));
numT = numel(trueExample.model.operators.T);
sizeT = size(trueExample.model.operators.T);

transmissibilityData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    transmissibilityData{i}.T = exp(logMean-0.5 + 2*logStd*randn(sizeT));
end

transSample = OperatorSamplesHM('data', transmissibilityData)

%% Select quantity of interest class matching the what we have as observations
% We validate the QoI with the trueProblem, since this will be our ensemble
% base problem as well.

qoi = WellQoIHM('wellNames', {'P1', 'P2'}, ...
              'fldname', {'qOs', 'qWs'}, ...
              'observationResultHandler', observationResultHandler, ...
              'observationCov', obsStdDev^2);


%% Create the ensemble
transEnsemble = MRSTHistoryMatchingEnsemble(trueExample, transSample, qoi, ...
    'directory', fullfile(topDirectory, 'trans'), ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true)

%% Displaying the observations and observation error cov through the ensemble
disp('observation and scaling vector')
[obsVector, obsScaling] = transEnsemble.qoi.getObservationAndScaling()
disp('observation error covariance matrix')
size(transEnsemble.qoi.getObservationErrorCov())

%% Run ensemble
transEnsemble.simulateEnsembleMembers();

%% Get simulated observations
disp('simulated observations')
size(transEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
size(transEnsemble.getEnsembleSamples())

%% Do history matching and thereby update the samples in the ensemble
disp('updated sample object:')
transEnsemble.doHistoryMatching()


%% Run ensemble with new samples
transEnsemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results

transEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'legend', {'observations', 'posterior mean', 'prior mean'});






%% Create another ensemble using stochastic porovolumes
% ---------------------

pvData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    poro = gaussianField(trueExample.model.G.cartDims, [0.1 0.7]); 
    pvData{i}.pv = poro(:).*trueExample.model.G.cells.volumes;
end

pvSamples = OperatorSamplesHM('data', pvData)



%% Define new ensemble
pvEnsemble = MRSTHistoryMatchingEnsemble(trueExample, pvSamples, qoi, ...
    'directory', fullfile(topDirectory, 'pv'), ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Simulate
pvEnsemble.simulateEnsembleMembers();

%% Do history matching
pvEnsemble.doHistoryMatching()

%% Run new ensemble
pvEnsemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results
pvEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'legend', {'observations', 'posterior mean', 'prior mean'});




%% Create another ensemble using stochastic porovolumes and trans
% ---------------------

operatorData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    operatorData{i}.pv = pvData{i}.pv;
    operatorData{i}.T  = transmissibilityData{i}.T;
end

operatorSamples = OperatorSamplesHM('data', operatorData)



%% Define new ensemble
operatorEnsemble = MRSTHistoryMatchingEnsemble(trueExample, operatorSamples, qoi, ...
    'directory', fullfile(topDirectory, 'operator'), ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Simulate
operatorEnsemble.simulateEnsembleMembers();

%% Do history matching
operatorEnsemble.doHistoryMatching()

%% Run new ensemble
operatorEnsemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results
operatorEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'legend', {'observations', 'posterior mean', 'prior mean'});











