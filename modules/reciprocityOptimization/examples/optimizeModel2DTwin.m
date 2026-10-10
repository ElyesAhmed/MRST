%% optimizeModel2D - NPV Optimization for 2D Reservoir Model
% Optimize well controls for maximum NPV using three different optimization strategies

mrstModule add ad-core ad-blackoil ad-props optimization spe10
setupModel2D

%% Initialize Reservoir Model
model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
state0 = initResSol(G, 200*barsa, [0, 1]);

%% NPV Objective Function Parameters
discountRate = 0.9;      % Yearly discount factor
oilPrice = 60;            % Oil revenue ($/stb)
waterProdCost = 6;        % Water production handling costs ($/stb)
waterInjCost = 6;         % Water injection cost ($/stb)

% Define two sets of NPV options for alternating optimization
npvOptionsDiscounted = {'OilPrice', oilPrice, ...
                       'WaterProductionCost', waterProdCost, ...
                       'WaterInjectionCost', waterInjCost, ...
                       'DiscountFactor', discountRate};

npvOptionsUndiscounted = {'OilPrice', oilPrice, ...
                         'WaterProductionCost', waterProdCost, ...
                         'WaterInjectionCost', waterInjCost, ...
                         'DiscountFactor', 0}; % No discounting

%% Optimization 1: RATE Controls for Injectors, BHP Controls for Producers
fprintf('=== Optimization 1: RATE-Injectors / BHP-Producers ===\n');

% Define control bounds for rate-controlled injectors and BHP-controlled producers
injRateLimits = [0 600]/day;          % Injector rate limits
prodBHPLimits = [10 500]*barsa;      % Producer BHP limits

scalingRateBhp = struct();
scalingRateBhp.boxLims = [injRateLimits; injRateLimits; prodBHPLimits; prodBHPLimits];
scalingRateBhp.obj = 3.2e7;           % Objective function scaling

% Get initial scaled controls
initialControlsRateBhp = schedule2control(schedule, scalingRateBhp);

% Define objective function with discounting
objectiveFunctionDiscounted = @(model, states, sched, varargin) ...
    NPVOW(model, states, sched, varargin{:}, npvOptionsDiscounted{:});

% Create evaluation function handle
evalNpvRateBhp = @(controls) evalObjective(controls, objectiveFunctionDiscounted, ...
    state0, model, schedule, scalingRateBhp);

% Run optimization
[optimalValueRateBhp, optimalControlsRateBhp, historyRateBhp] = ...
    unitBoxBFGS(initialControlsRateBhp, evalNpvRateBhp, 'maxIt', 25);

% Create optimized schedule
optimizedScheduleRateBhp = control2schedule(optimalControlsRateBhp, schedule, scalingRateBhp);

% Save optimized schedule
outputPath = fullfile(mrstPath('optimization'), 'examples', 'model2D', 'optimized_schedule_rate_bhp.mat');
save(outputPath, 'optimizedScheduleRateBhp');

% Simulate optimized schedule
[wellSolutionsRateBhp, statesRateBhp] = simulateScheduleAD(state0, model, optimizedScheduleRateBhp);
npvRateBhp = NPVOW(model, statesRateBhp, optimizedScheduleRateBhp, npvOptionsDiscounted{:});

%% Convert to Alternative Control Strategy: BHP Injectors / RATE Producers
fprintf('=== Converting to BHP-Injectors / RATE-Producers Strategy ===\n');

controlConversionMap = struct('bhp', 'rate', 'rate', 'bhp');
scheduleBhpRate = convertWellControls(optimizedScheduleRateBhp, statesRateBhp, model, ...
    'ConversionMap', controlConversionMap, 'bc', []);

%% Optimization 2: BHP Controls for Injectors, RATE Controls for Producers
fprintf('=== Optimization 2: BHP-Injectors / RATE-Producers ===\n');

% Define control bounds for BHP-controlled injectors and rate-controlled producers
injBHPLimits = [20 3000]*barsa;      % Injector BHP limits
prodRateLimits = [-600 0]/day;        % Producer rate limits (negative for production)

scalingBhpRate = struct();
scalingBhpRate.boxLims = [injBHPLimits; injBHPLimits; prodRateLimits; prodRateLimits];
scalingBhpRate.obj = 1;               % Objective function scaling

% Get initial scaled controls for new strategy
initialControlsBhpRate = schedule2control(scheduleBhpRate, scalingBhpRate);

% Create evaluation function for new control strategy
evalNpvBhpRate = @(controls) evalObjective(controls, objectiveFunctionDiscounted, ...
    state0, model, scheduleBhpRate, scalingBhpRate);

% Run optimization for new control strategy
[optimalValueBhpRate, optimalControlsBhpRate, historyBhpRate] = ...
    unitBoxBFGS(initialControlsBhpRate, evalNpvBhpRate, 'maxIt', 25);

% Create optimized schedule for new strategy
optimizedScheduleBhpRate = control2schedule(optimalControlsBhpRate, scheduleBhpRate, scalingBhpRate);

% Simulate optimized schedule for new strategy
[wellSolutionsBhpRate, statesBhpRate] = simulateScheduleAD(state0, model, optimizedScheduleBhpRate);
npvBhpRate = NPVOW(model, statesBhpRate, optimizedScheduleBhpRate, npvOptionsDiscounted{:});

%% Optimization 3: Alternating Discount Factor Optimization
fprintf('=== Optimization 3: Alternating Discount Factor Strategy ===\n');

% Initialize for alternating optimization
% initialControlsBhpRate = initialControlsBhpRate; % Start from previous optimization
maxAlternations = 5;                     % Number of alternations
mixingParameter = 0.7;                   % Control mixing parameter

alternationHistory = struct();
alternationHistory.discountedNpv = [];
alternationHistory.undiscountedNpv = [];
alternationHistory.controls = [];

fprintf('Starting alternating optimization with %d alternations...\n', maxAlternations);

for altIter = 1:maxAlternations
    fprintf('\n--- Alternation %d/%d ---\n', altIter, maxAlternations);
    
    % Phase 1: Optimize with discounting (b = 0.05)
    fprintf('Phase 1: Optimizing with discount factor 0.05\n');
    evalDiscounted = @(controls) evalObjective(controls, objectiveFunctionDiscounted, ...
        state0, model, scheduleBhpRate, scalingBhpRate);
    
    [~, discountedControls, historyDiscounted] = unitBoxBFGS(initialControlsBhpRate, evalDiscounted, 'maxIt', 10);
        
    discountedControls = mixingParameter * discountedControls + (1 - mixingParameter) * initialControlsBhpRate;

    optimizedScheduleBhpRate = control2schedule(discountedControls, scheduleBhpRate, scalingBhpRate);
    [wellSolutionsBhpRate, statesBhpRate] = simulateScheduleAD(state0, model, optimizedScheduleBhpRate);

    scheduleRateBhp = convertWellControls(optimizedScheduleBhpRate, statesBhpRate, model, ...
    'ConversionMap', controlConversionMap, 'bc', []);
    
    initialControlsRateBhp = schedule2control(scheduleRateBhp, scalingRateBhp);

    % Phase 2: Optimize without discounting (b = 0)
    fprintf('Phase 2: Optimizing with discount factor 0\n');
    objectiveFunctionUndiscounted = @(model, states, sched, varargin) ...
        NPVOW(model, states, sched, varargin{:}, npvOptionsUndiscounted{:});
    
    evalUndiscounted = @(controls) evalObjective(controls, objectiveFunctionUndiscounted, ...
        state0, model, scheduleRateBhp, scalingRateBhp);
    
    [~, undiscountedControls, historyUndiscounted] = unitBoxBFGS(initialControlsRateBhp, evalUndiscounted, 'maxIt', 10);


    undiscountedControls = mixingParameter * undiscountedControls + (1 - mixingParameter) * initialControlsRateBhp;


      optimizedScheduleRateBhp = control2schedule(undiscountedControls, scheduleRateBhp, scalingRateBhp);
    [wellSolutionsRateBhp, statesRateBhp] = simulateScheduleAD(state0, model, optimizedScheduleRateBhp);

    scheduleBhpRate = convertWellControls(optimizedScheduleRateBhp, statesRateBhp, model, ...
    'ConversionMap', controlConversionMap, 'bc', []);

    initialControlsBhpRate = schedule2control(scheduleBhpRate, scalingBhpRate);

        
    
    % Evaluate both objectives with mixed controls
    mixedNpvDiscounted = evalDiscounted(initialControlsBhpRate);
    mixedNpvUndiscounted = evalUndiscounted(initialControlsBhpRate);
    
    % Store history
    alternationHistory.discountedNpv(end+1) = mixedNpvDiscounted;
    alternationHistory.undiscountedNpv(end+1) = mixedNpvUndiscounted;
    alternationHistory.controls(:, end+1) = initialControlsBhpRate;
    
    fprintf('Mixed NPV (discounted): $%.2f\n', mixedNpvDiscounted);
    fprintf('Mixed NPV (undiscounted): $%.2f\n', mixedNpvUndiscounted);
end

% Create final optimized schedule from alternating optimization
optimizedScheduleAlternating = control2schedule(initialControlsBhpRate, scheduleBhpRate, scalingBhpRate);

% Simulate final optimized schedule
[wellSolutionsAlternating, statesAlternating] = simulateScheduleAD(state0, model, optimizedScheduleAlternating);
npvAlternatingDiscounted = NPVOW(model, statesAlternating, optimizedScheduleAlternating, npvOptionsDiscounted{:});
npvAlternatingUndiscounted = NPVOW(model, statesAlternating, optimizedScheduleAlternating, npvOptionsUndiscounted{:});

%% Results Summary
fprintf('\n=== Optimization Results Summary ===\n');
fprintf('Strategy 1 (RATE-Inj / BHP-Prod): NPV = $%.2f\n', sum(vertcat(npvRateBhp{:})));
fprintf('Strategy 2 (BHP-Inj / RATE-Prod): NPV = $%.2f\n', sum(vertcat(npvBhpRate{:})));
fprintf('Strategy 3 (Alternating DF): NPV = $%.2f (discounted), $%.2f (undiscounted)\n', ...
    sum(vertcat(npvAlternatingDiscounted{:})), sum(vertcat(npvAlternatingUndiscounted{:})));

fprintf('\nImprovement over Strategy 1: %.2f%%\n', ...
    (sum(vertcat(npvAlternatingDiscounted{:}))/sum(vertcat(npvRateBhp{:})) - 1)*100);

%% Plot Optimization History
figure;

% Strategy 1 & 2 NPV progression
subplot(3,1,1);
timeSteps = 1:numel(vertcat(npvRateBhp{:}));
plot(timeSteps, cumsum(vertcat(npvRateBhp{:})), 'b-o', 'LineWidth', 2, 'DisplayName', 'RATE-Inj/BHP-Prod');
hold on;
plot(timeSteps, cumsum(vertcat(npvBhpRate{:})), 'r-s', 'LineWidth', 2, 'DisplayName', 'BHP-Inj/RATE-Prod');
plot(timeSteps, cumsum(vertcat(npvAlternatingDiscounted{:})), 'g-^', 'LineWidth', 2, 'DisplayName', 'Alternating DF');
title('Cumulative NPV Comparison');
xlabel('Time Step');
ylabel('Cumulative NPV ($)');
legend('show');
grid on;

% Alternation history
subplot(3,1,2);
plot(1:maxAlternations, alternationHistory.discountedNpv, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
plot(1:maxAlternations, alternationHistory.undiscountedNpv, 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
title('Alternating Optimization Progress');
xlabel('Alternation Number');
ylabel('NPV ($)');
legend('Discounted (b=0.05)', 'Undiscounted (b=0)');
grid on;

% Control evolution
subplot(3,1,3);
plot(alternationHistory.controls', 'LineWidth', 1.5);
title('Control Parameter Evolution During Alternating Optimization');
xlabel('Alternation Number');
ylabel('Control Value');
grid on;

%% Save Final Results
finalResults = struct();
finalResults.optimizedScheduleRateBhp = optimizedScheduleRateBhp;
finalResults.optimizedScheduleBhpRate = optimizedScheduleBhpRate;
finalResults.optimizedScheduleAlternating = optimizedScheduleAlternating;
finalResults.npvRateBhp = sum(vertcat(npvRateBhp{:}));
finalResults.npvBhpRate = sum(vertcat(npvBhpRate{:}));
finalResults.npvAlternatingDiscounted = sum(vertcat(npvAlternatingDiscounted{:}));
finalResults.npvAlternatingUndiscounted = sum(vertcat(npvAlternatingUndiscounted{:}));
finalResults.historyRateBhp = historyRateBhp;
finalResults.historyBhpRate = historyBhpRate;
finalResults.alternationHistory = alternationHistory;

save('final_optimization_results.mat', 'finalResults');

fprintf('\nAll optimizations completed successfully. Results saved.\n');
fprintf('Final discounted NPV: $%.2f\n', finalResults.npvAlternatingDiscounted);
fprintf('Final undiscounted NPV: $%.2f\n', finalResults.npvAlternatingUndiscounted);

%% Copyright Notice
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with MRST. If not, see <http://www.gnu.org/licenses/>.