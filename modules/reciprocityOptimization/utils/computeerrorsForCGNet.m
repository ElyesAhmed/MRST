%% Compute percentage errors for BHP, qOs, and qWs using proper time weighting

% Function to extract well data from cell array of well solutions
getWellData = @(ws) struct(...
    'bhp', cell2mat(cellfun(@(x) [x.bhp], ws, 'UniformOutput', false)), ...
    'qOs', cell2mat(cellfun(@(x) [x.qOs], ws, 'UniformOutput', false)), ...
    'qWs', cell2mat(cellfun(@(x) [x.qWs], ws, 'UniformOutput', false)) ...
);

% Get time steps from schedule
dt = schedule.step.val;
total_time = sum(dt);

% Define thresholds for zero protection (adjust as needed)
bhp_threshold = 1 * barsa;      % 1 bar threshold for BHP
rate_threshold = 1e-6 * (1/day); % Small rate threshold for qOs and qWs

% Function to calculate time-weighted percentage error for individual wells with zero protection
calculateWellError = @(opt_vals, ref_vals, dt, threshold) ...
    100 * sum(abs(opt_vals - ref_vals) .* dt, 1) ./ ...
    max(sum(abs(ref_vals) .* dt, 1), threshold * total_time);


calculateNormWellError = @(opt_vals, ref_vals, dt, threshold) ...
    100 * sum(abs(opt_vals - ref_vals) .* dt, 1) ./ ...
    max(range(ref_vals, 1) * total_time, threshold * total_time);

% Function to calculate total error (all wells combined) with zero protection
calculateTotalError = @(opt_vals, ref_vals, dt, threshold) ...
    100 * sum(abs(opt_vals - ref_vals) .* dt, 'all') ./ ...
    max(sum(abs(ref_vals) .* dt, 'all'), threshold * total_time * size(ref_vals, 2));
% Function to calculate total error (all wells combined) with zero protection
calculateNormTotalError = @(opt_vals, ref_vals, dt, threshold) ...
    100 * sum(abs(opt_vals - ref_vals) .* dt, 'all') ./ ...
    max(range(ref_vals, 'all') * total_time * size(ref_vals, 2), threshold);


% Extract reference and optimized well data
if NashAlgorithm
    % For Nash optimization case
    ref_data1 = getWellData(wsRef);
    opt_data1 = getWellData(wellSols_opt1);
    ref_data2 = getWellData(wsRef_new); 
    opt_data2 = getWellData(wellSols_opt2);
    
    % Calculate percentage errors for RATE-BHP problem
    bhp_error1 = calculateNormWellError(opt_data1.bhp, ref_data1.bhp, dt, bhp_threshold);
    qOs_error1 = calculateWellError(opt_data1.qOs, ref_data1.qOs, dt, rate_threshold);
    qWs_error1 = calculateWellError(opt_data1.qWs, ref_data1.qWs, dt, rate_threshold);
    
    % Calculate percentage errors for BHP-RATE problem
    bhp_error2 = calculateNormWellError(opt_data2.bhp, ref_data2.bhp, dt, bhp_threshold);
    qOs_error2 = calculateWellError(opt_data2.qOs, ref_data2.qOs, dt, rate_threshold);
    qWs_error2 = calculateWellError(opt_data2.qWs, ref_data2.qWs, dt, rate_threshold);
    
    % Calculate total errors (combining all wells)
    total_bhp_error1 = calculateNormTotalError(opt_data1.bhp, ref_data1.bhp, dt, bhp_threshold);
    total_qOs_error1 = calculateTotalError(opt_data1.qOs, ref_data1.qOs, dt, rate_threshold);
    total_qWs_error1 = calculateTotalError(opt_data1.qWs, ref_data1.qWs, dt, rate_threshold);
    
    total_bhp_error2 = calculateNormTotalError(opt_data2.bhp, ref_data2.bhp, dt, bhp_threshold);
    total_qOs_error2 = calculateTotalError(opt_data2.qOs, ref_data2.qOs, dt, rate_threshold);
    total_qWs_error2 = calculateTotalError(opt_data2.qWs, ref_data2.qWs, dt, rate_threshold);
    
    % Handle cases where reference values are zero (will produce Inf/NaN)
    bhp_error1(isinf(bhp_error1) | isnan(bhp_error1)) = 0;
    qOs_error1(isinf(qOs_error1) | isnan(qOs_error1)) = 0;
    qWs_error1(isinf(qWs_error1) | isnan(qWs_error1)) = 0;
    
    bhp_error2(isinf(bhp_error2) | isnan(bhp_error2)) = 0;
    qOs_error2(isinf(qOs_error2) | isnan(qOs_error2)) = 0;
    qWs_error2(isinf(qWs_error2) | isnan(qWs_error2)) = 0;
    
    % Display results for RATE-BHP problem
    fprintf('\n=== RATE-BHP Problem - Time-Weighted Percentage Errors ===\n');
    fprintf('Total simulation time: %.2f days\n', total_time/day);
    fprintf('Rate threshold: %.2e m³/day\n', rate_threshold*(1/day));
    fprintf('BHP Mean Error: %.2f%%, Max Error: %.2f%%\n', mean(bhp_error1, 'omitnan'), max(bhp_error1, [], 'omitnan'));
    fprintf('qOs Mean Error: %.2f%%, Max Error: %.2f%%\n', mean(qOs_error1, 'omitnan'), max(qOs_error1, [], 'omitnan'));
    fprintf('qWs Mean Error: %.2f%%, Max Error: %.2f%%\n', mean(qWs_error1, 'omitnan'), max(qWs_error1, [], 'omitnan'));
    fprintf('\n--- TOTAL ERRORS (All Wells Combined) ---\n');
    fprintf('Total BHP Error: %.2f%%\n', total_bhp_error1);
    fprintf('Total qOs Error: %.2f%%\n', total_qOs_error1);
    fprintf('Total qWs Error: %.2f%%\n', total_qWs_error1);
    
    % Display results for BHP-RATE problem
    fprintf('\n=== BHP-RATE Problem - Time-Weighted Percentage Errors ===\n');
    fprintf('BHP Mean Error: %.2f%%, Max Error: %.2f%%\n', mean(bhp_error2, 'omitnan'), max(bhp_error2, [], 'omitnan'));
    fprintf('qOs Mean Error: %.2f%%, Max Error: %.2f%%\n', mean(qOs_error2, 'omitnan'), max(qOs_error2, [], 'omitnan'));
    fprintf('qWs Mean Error: %.2f%%, Max Error: %.2f%%\n', mean(qWs_error2, 'omitnan'), max(qWs_error2, [], 'omitnan'));
    fprintf('\n--- TOTAL ERRORS (All Wells Combined) ---\n');
    fprintf('Total BHP Error: %.2f%%\n', total_bhp_error2);
    fprintf('Total qOs Error: %.2f%%\n', total_qOs_error2);
    fprintf('Total qWs Error: %.2f%%\n', total_qWs_error2);
    
else
    % For single optimization case
    ref_data = getWellData(wsRef);
    if playerChoice == 1
        opt_data = getWellData(wellSols_opt1);
    else
        opt_data = getWellData(wellSols_opt2);
    end

    % Calculate percentage errors with zero protection
    bhp_error = calculateNormWellError(opt_data.bhp, ref_data.bhp, dt, bhp_threshold);
    qOs_error = calculateWellError(opt_data.qOs, ref_data.qOs, dt, rate_threshold);
    qWs_error = calculateWellError(opt_data.qWs, ref_data.qWs, dt, rate_threshold);
    
    % Calculate total errors (combining all wells)
    total_bhp_error = calculateNormTotalError(opt_data.bhp, ref_data.bhp, dt, bhp_threshold);
    total_qOs_error = calculateTotalError(opt_data.qOs, ref_data.qOs, dt, rate_threshold);
    total_qWs_error = calculateTotalError(opt_data.qWs, ref_data.qWs, dt, rate_threshold);
    
    % Handle cases where reference values are zero (will produce Inf/NaN)
    bhp_error(isinf(bhp_error) | isnan(bhp_error)) = 0;
    qOs_error(isinf(qOs_error) | isnan(qOs_error)) = 0;
    qWs_error(isinf(qWs_error) | isnan(qWs_error)) = 0;
    
    % Display results
    fprintf('\n=== RATE-BHP Problem - Time-Weighted Percentage Errors ===\n');
    fprintf('Total simulation time: %.2f days\n', total_time/day);
    fprintf('Rate threshold: %.2e m³/day\n', rate_threshold*(1/day));
    fprintf('BHP Mean Error: %.2f%%, Max Error: %.2f%%\n', mean(bhp_error, 'omitnan'), max(bhp_error, [], 'omitnan'));
    fprintf('qOs Mean Error: %.2f%%, Max Error: %.2f%%\n', mean(qOs_error, 'omitnan'), max(qOs_error, [], 'omitnan'));
    fprintf('qWs Mean Error: %.2f%%, Max Error: %.2f%%\n', mean(qWs_error, 'omitnan'), max(qWs_error, [], 'omitnan'));
    fprintf('\n--- TOTAL ERRORS (All Wells Combined) ---\n');
    fprintf('Total BHP Error: %.2f%%\n', total_bhp_error);
    fprintf('Total qOs Error: %.2f%%\n', total_qOs_error);
    fprintf('Total qWs Error: %.2f%%\n', total_qWs_error);
end

%% Additional analysis: Identify wells with zero or very small reference values
if NashAlgorithm
    fprintf('\n=== Analysis of Wells with Zero/Small Reference Values (RATE-BHP) ===\n');
    num_wells = size(ref_data1.bhp, 2);
    small_well_count = 0;
    
    for well_idx = 1:num_wells
        total_ref_bhp = sum(abs(ref_data1.bhp(:, well_idx)) .* dt);
        total_ref_qOs = sum(abs(ref_data1.qOs(:, well_idx)) .* dt);
        total_ref_qWs = sum(abs(ref_data1.qWs(:, well_idx)) .* dt);
        
        if total_ref_qOs < rate_threshold * total_time || total_ref_qWs < rate_threshold * total_time
            small_well_count = small_well_count + 1;
            fprintf('Well %d has small integrated reference rates:\n', well_idx);
            fprintf('  Integrated qOs reference: %.2e m³ (threshold: %.2e m³)\n', total_ref_qOs, rate_threshold * total_time);
            fprintf('  Integrated qWs reference: %.2e m³ (threshold: %.2e m³)\n', total_ref_qWs, rate_threshold * total_time);
            fprintf('  qOs error (protected): %.2f%%\n', qOs_error1(well_idx));
            fprintf('  qWs error (protected): %.2f%%\n', qWs_error1(well_idx));
        end
    end
    fprintf('Total wells with small reference rates: %d/%d\n', small_well_count, num_wells);
end

%% Well-by-well error breakdown
if NashAlgorithm
    fprintf('\n=== Well-by-Well Error Analysis (RATE-BHP) ===\n');
    num_wells = size(ref_data1.bhp, 2);
    for well_idx = 1:num_wells
        fprintf('Well %d:\n', well_idx);
        fprintf('  BHP error: %.2f%%\n', bhp_error1(well_idx));
        fprintf('  qOs error: %.2f%%\n', qOs_error1(well_idx));
        fprintf('  qWs error: %.2f%%\n', qWs_error1(well_idx));
        
        % Show time-integrated reference values for context
        int_bhp = sum(abs(ref_data1.bhp(:, well_idx)) .* dt)/barsa;
        int_qOs = sum(abs(ref_data1.qOs(:, well_idx)) .* dt);
        int_qWs = sum(abs(ref_data1.qWs(:, well_idx)) .* dt);
        fprintf('  (Integrated reference: BHP=%.1f bar·day, qOs=%.2e m³, qWs=%.2e m³)\n', ...
                int_bhp, int_qOs, int_qWs);
    end
end

%% Compute validation errors with zero protection
fprintf('\n\n=== VALIDATION SCHEDULE RESULTS ===\n');

if NashAlgorithm
    %% Nash Algorithm Case - Calculate errors for both players
    
    % Player 1 (RATE-BHP) validation errors
    val_ref_data1 = getWellData(ws_new1);
    val_coarse_old_data1 = getWellData(ws_coarse_old1);
    val_coarse_opt_data1 = getWellData(ws_coarse_new1);
    
    % Player 2 (BHP-RATE) validation errors  
    val_ref_data2 = getWellData(ws_new2);
    val_coarse_old_data2 = getWellData(ws_coarse_old2);
    val_coarse_opt_data2 = getWellData(ws_coarse_new2);
    
    % Calculate errors for Player 1 with zero protection
    val_bhp_error_old1 = calculateNormWellError(val_coarse_old_data1.bhp, val_ref_data1.bhp, dt, bhp_threshold);
    val_qOs_error_old1 = calculateWellError(val_coarse_old_data1.qOs, val_ref_data1.qOs, dt, rate_threshold);
    val_qWs_error_old1 = calculateWellError(val_coarse_old_data1.qWs, val_ref_data1.qWs, dt, rate_threshold);
    
    val_bhp_error_opt1 = calculateNormWellError(val_coarse_opt_data1.bhp, val_ref_data1.bhp, dt, bhp_threshold);
    val_qOs_error_opt1 = calculateWellError(val_coarse_opt_data1.qOs, val_ref_data1.qOs, dt, rate_threshold);
    val_qWs_error_opt1 = calculateWellError(val_coarse_opt_data1.qWs, val_ref_data1.qWs, dt, rate_threshold);
    
    val_total_bhp_error_old1 = calculateNormTotalError(val_coarse_old_data1.bhp, val_ref_data1.bhp, dt, bhp_threshold);
    val_total_qOs_error_old1 = calculateTotalError(val_coarse_old_data1.qOs, val_ref_data1.qOs, dt, rate_threshold);
    val_total_qWs_error_old1 = calculateTotalError(val_coarse_old_data1.qWs, val_ref_data1.qWs, dt, rate_threshold);
    
    val_total_bhp_error_opt1 = calculateNormTotalError(val_coarse_opt_data1.bhp, val_ref_data1.bhp, dt, bhp_threshold);
    val_total_qOs_error_opt1 = calculateTotalError(val_coarse_opt_data1.qOs, val_ref_data1.qOs, dt, rate_threshold);
    val_total_qWs_error_opt1 = calculateTotalError(val_coarse_opt_data1.qWs, val_ref_data1.qWs, dt, rate_threshold);
    
    % Calculate errors for Player 2 with zero protection
    val_bhp_error_old2 = calculateNormWellError(val_coarse_old_data2.bhp, val_ref_data2.bhp, dt, bhp_threshold);
    val_qOs_error_old2 = calculateWellError(val_coarse_old_data2.qOs, val_ref_data2.qOs, dt, rate_threshold);
    val_qWs_error_old2 = calculateWellError(val_coarse_old_data2.qWs, val_ref_data2.qWs, dt, rate_threshold);
    
    val_bhp_error_opt2 = calculateNormWellError(val_coarse_opt_data2.bhp, val_ref_data2.bhp, dt, bhp_threshold);
    val_qOs_error_opt2 = calculateWellError(val_coarse_opt_data2.qOs, val_ref_data2.qOs, dt, rate_threshold);
    val_qWs_error_opt2 = calculateWellError(val_coarse_opt_data2.qWs, val_ref_data2.qWs, dt, rate_threshold);
    
    val_total_bhp_error_old2 = calculateNormTotalError(val_coarse_old_data2.bhp, val_ref_data2.bhp, dt, bhp_threshold);
    val_total_qOs_error_old2 = calculateTotalError(val_coarse_old_data2.qOs, val_ref_data2.qOs, dt, rate_threshold);
    val_total_qWs_error_old2 = calculateTotalError(val_coarse_old_data2.qWs, val_ref_data2.qWs, dt, rate_threshold);
    
    val_total_bhp_error_opt2 = calculateNormTotalError(val_coarse_opt_data2.bhp, val_ref_data2.bhp, dt, bhp_threshold);
    val_total_qOs_error_opt2 = calculateTotalError(val_coarse_opt_data2.qOs, val_ref_data2.qOs, dt, rate_threshold);
    val_total_qWs_error_opt2 = calculateTotalError(val_coarse_opt_data2.qWs, val_ref_data2.qWs, dt, rate_threshold);
    
    % Handle cases where reference values are zero for both players
    % Player 1
    val_bhp_error_old1(isinf(val_bhp_error_old1) | isnan(val_bhp_error_old1)) = 0;
    val_qOs_error_old1(isinf(val_qOs_error_old1) | isnan(val_qOs_error_old1)) = 0;
    val_qWs_error_old1(isinf(val_qWs_error_old1) | isnan(val_qWs_error_old1)) = 0;
    val_bhp_error_opt1(isinf(val_bhp_error_opt1) | isnan(val_bhp_error_opt1)) = 0;
    val_qOs_error_opt1(isinf(val_qOs_error_opt1) | isnan(val_qOs_error_opt1)) = 0;
    val_qWs_error_opt1(isinf(val_qWs_error_opt1) | isnan(val_qWs_error_opt1)) = 0;
    
    % Player 2
    val_bhp_error_old2(isinf(val_bhp_error_old2) | isnan(val_bhp_error_old2)) = 0;
    val_qOs_error_old2(isinf(val_qOs_error_old2) | isnan(val_qOs_error_old2)) = 0;
    val_qWs_error_old2(isinf(val_qWs_error_old2) | isnan(val_qWs_error_old2)) = 0;
    val_bhp_error_opt2(isinf(val_bhp_error_opt2) | isnan(val_bhp_error_opt2)) = 0;
    val_qOs_error_opt2(isinf(val_qOs_error_opt2) | isnan(val_qOs_error_opt2)) = 0;
    val_qWs_error_opt2(isinf(val_qWs_error_opt2) | isnan(val_qWs_error_opt2)) = 0;
    
    % Display results for both players
    fprintf('\n=== VALIDATION - Player 1 (RATE-BHP) ===\n');
    fprintf('Initial Model - Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
            val_total_bhp_error_old1, val_total_qOs_error_old1, val_total_qWs_error_old1);
    fprintf('Optimized Model - Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
            val_total_bhp_error_opt1, val_total_qOs_error_opt1, val_total_qWs_error_opt1);
    
    fprintf('\n=== VALIDATION - Player 2 (BHP-RATE) ===\n');
    fprintf('Initial Model - Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
            val_total_bhp_error_old2, val_total_qOs_error_old2, val_total_qWs_error_old2);
    fprintf('Optimized Model - Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
            val_total_bhp_error_opt2, val_total_qOs_error_opt2, val_total_qWs_error_opt2);
    
    % Calculate improvement percentages for both players
    bhp_improvement1 = ((val_total_bhp_error_old1 - val_total_bhp_error_opt1) / max(val_total_bhp_error_old1, 1e-10)) * 100;
    qOs_improvement1 = ((val_total_qOs_error_old1 - val_total_qOs_error_opt1) / max(val_total_qOs_error_old1, 1e-10)) * 100;
    qWs_improvement1 = ((val_total_qWs_error_old1 - val_total_qWs_error_opt1) / max(val_total_qWs_error_old1, 1e-10)) * 100;
    
    bhp_improvement2 = ((val_total_bhp_error_old2 - val_total_bhp_error_opt2) / max(val_total_bhp_error_old2, 1e-10)) * 100;
    qOs_improvement2 = ((val_total_qOs_error_old2 - val_total_qOs_error_opt2) / max(val_total_qOs_error_old2, 1e-10)) * 100;
    qWs_improvement2 = ((val_total_qWs_error_old2 - val_total_qWs_error_opt2) / max(val_total_qWs_error_old2, 1e-10)) * 100;
    
    fprintf('\n=== VALIDATION - Improvement Over Initial ===\n');
    fprintf('Player 1 - BHP: %.1f%%, qOs: %.1f%%, qWs: %.1f%% reduction\n', ...
            max(0, bhp_improvement1), max(0, qOs_improvement1), max(0, qWs_improvement1));
    fprintf('Player 2 - BHP: %.1f%%, qOs: %.1f%%, qWs: %.1f%% reduction\n', ...
            max(0, bhp_improvement2), max(0, qOs_improvement2), max(0, qWs_improvement2));
    
else
    %% Non-Nash Case - Calculate errors for the chosen player only
    
    if playerChoice == 1
        % Player 1 (RATE-BHP) validation
        val_ref_data = getWellData(ws_new1);
        val_coarse_old_data = getWellData(ws_coarse_old1);
        val_coarse_opt_data = getWellData(ws_coarse_new1);
        player_name = 'Player 1 (RATE-BHP)';
    else
        % Player 2 (BHP-RATE) validation
        val_ref_data = getWellData(ws_new2);
        val_coarse_old_data = getWellData(ws_coarse_old2);
        val_coarse_opt_data = getWellData(ws_coarse_new2);
        player_name = 'Player 2 (BHP-RATE)';
    end
    
    % Calculate percentage errors with zero protection
    val_bhp_error_old = calculateNormWellError(val_coarse_old_data.bhp, val_ref_data.bhp, dt, bhp_threshold);
    val_qOs_error_old = calculateWellError(val_coarse_old_data.qOs, val_ref_data.qOs, dt, rate_threshold);
    val_qWs_error_old = calculateWellError(val_coarse_old_data.qWs, val_ref_data.qWs, dt, rate_threshold);
    
    val_bhp_error_opt = calculateNormWellError(val_coarse_opt_data.bhp, val_ref_data.bhp, dt, bhp_threshold);
    val_qOs_error_opt = calculateWellError(val_coarse_opt_data.qOs, val_ref_data.qOs, dt, rate_threshold);
    val_qWs_error_opt = calculateWellError(val_coarse_opt_data.qWs, val_ref_data.qWs, dt, rate_threshold);
    
    % Calculate total errors
    val_total_bhp_error_old = calculateNormTotalError(val_coarse_old_data.bhp, val_ref_data.bhp, dt, bhp_threshold);
    val_total_qOs_error_old = calculateTotalError(val_coarse_old_data.qOs, val_ref_data.qOs, dt, rate_threshold);
    val_total_qWs_error_old = calculateTotalError(val_coarse_old_data.qWs, val_ref_data.qWs, dt, rate_threshold);
    
    val_total_bhp_error_opt = calculateNormTotalError(val_coarse_opt_data.bhp, val_ref_data.bhp, dt, bhp_threshold);
    val_total_qOs_error_opt = calculateTotalError(val_coarse_opt_data.qOs, val_ref_data.qOs, dt, rate_threshold);
    val_total_qWs_error_opt = calculateTotalError(val_coarse_opt_data.qWs, val_ref_data.qWs, dt, rate_threshold);
    
    % Handle cases where reference values are zero
    val_bhp_error_old(isinf(val_bhp_error_old) | isnan(val_bhp_error_old)) = 0;
    val_qOs_error_old(isinf(val_qOs_error_old) | isnan(val_qOs_error_old)) = 0;
    val_qWs_error_old(isinf(val_qWs_error_old) | isnan(val_qWs_error_old)) = 0;
    
    val_bhp_error_opt(isinf(val_bhp_error_opt) | isnan(val_bhp_error_opt)) = 0;
    val_qOs_error_opt(isinf(val_qOs_error_opt) | isnan(val_qOs_error_opt)) = 0;
    val_qWs_error_opt(isinf(val_qWs_error_opt) | isnan(val_qWs_error_opt)) = 0;
    
    % Display results
    fprintf('\n=== VALIDATION - %s ===\n', player_name);
    fprintf('Initial Model - Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
            val_total_bhp_error_old, val_total_qOs_error_old, val_total_qWs_error_old);
    fprintf('Optimized Model - Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
            val_total_bhp_error_opt, val_total_qOs_error_opt, val_total_qWs_error_opt);
    
    % Calculate improvement percentages
    bhp_improvement = ((val_total_bhp_error_old - val_total_bhp_error_opt) / max(val_total_bhp_error_old, 1e-10)) * 100;
    qOs_improvement = ((val_total_qOs_error_old - val_total_qOs_error_opt) / max(val_total_qOs_error_old, 1e-10)) * 100;
    qWs_improvement = ((val_total_qWs_error_old - val_total_qWs_error_opt) / max(val_total_qWs_error_old, 1e-10)) * 100;
    
    fprintf('\n=== VALIDATION - Improvement Over Initial ===\n');
    fprintf('BHP Improvement: %.1f%% reduction\n', max(0, bhp_improvement));
    fprintf('qOs Improvement: %.1f%% reduction\n', max(0, qOs_improvement));
    fprintf('qWs Improvement: %.1f%% reduction\n', max(0, qWs_improvement));
end

%% Well-by-well detailed error analysis table
fprintf('\n\n=== DETAILED WELL-BY-WELL ERROR ANALYSIS ===\n');

if NashAlgorithm
    %% Nash Algorithm - Show both models
    
    fprintf('\n--- MODEL 1 (RATE-BHP) - Well-by-Well Training Errors ---\n');
    fprintf('%-6s %-10s %-10s %-10s %-15s %-15s %-15s\n', ...
            'Well', 'BHP(%)', 'qOs(%)', 'qWs(%)', 'BHP_ref', 'qOs_ref', 'qWs_ref');
    fprintf('%-6s %-10s %-10s %-10s %-15s %-15s %-15s\n', ...
            '----', '------', '------', '------', '------', '------', '------');
    
    num_wells = size(ref_data1.bhp, 2);
    for well_idx = 1:num_wells
        int_bhp = sum(abs(ref_data1.bhp(:, well_idx)) .* dt)/barsa;
        int_qOs = sum(abs(ref_data1.qOs(:, well_idx)) .* dt);
        int_qWs = sum(abs(ref_data1.qWs(:, well_idx)) .* dt);
        
        fprintf('%-6d %-10.2f %-10.2f %-10.2f %-15.1f %-15.2e %-15.2e\n', ...
                well_idx, bhp_error1(well_idx), qOs_error1(well_idx), qWs_error1(well_idx), ...
                int_bhp, int_qOs, int_qWs);
    end
    
    fprintf('\n--- MODEL 2 (BHP-RATE) - Well-by-Well Training Errors ---\n');
    fprintf('%-6s %-10s %-10s %-10s %-15s %-15s %-15s\n', ...
            'Well', 'BHP(%)', 'qOs(%)', 'qWs(%)', 'BHP_ref', 'qOs_ref', 'qWs_ref');
    fprintf('%-6s %-10s %-10s %-10s %-15s %-15s %-15s\n', ...
            '----', '------', '------', '------', '------', '------', '------');
    
    for well_idx = 1:num_wells
        int_bhp = sum(abs(ref_data2.bhp(:, well_idx)) .* dt)/barsa;
        int_qOs = sum(abs(ref_data2.qOs(:, well_idx)) .* dt);
        int_qWs = sum(abs(ref_data2.qWs(:, well_idx)) .* dt);
        
        fprintf('%-6d %-10.2f %-10.2f %-10.2f %-15.1f %-15.2e %-15.2e\n', ...
                well_idx, bhp_error2(well_idx), qOs_error2(well_idx), qWs_error2(well_idx), ...
                int_bhp, int_qOs, int_qWs);
    end
    
    % Validation errors table
    fprintf('\n--- MODEL 1 (RATE-BHP) - Well-by-Well Validation Errors ---\n');
    fprintf('%-6s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
            'Well', 'BHP_Init(%)', 'BHP_Opt(%)', 'qOs_Init(%)', 'qOs_Opt(%)', 'qWs_Init(%)', 'qWs_Opt(%)');
    fprintf('%-6s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
            '----', '-----------', '-----------', '-----------', '-----------', '-----------', '-----------');
    
    for well_idx = 1:num_wells
        fprintf('%-6d %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f\n', ...
                well_idx, val_bhp_error_old1(well_idx), val_bhp_error_opt1(well_idx), ...
                val_qOs_error_old1(well_idx), val_qOs_error_opt1(well_idx), ...
                val_qWs_error_old1(well_idx), val_qWs_error_opt1(well_idx));
    end
    
    fprintf('\n--- MODEL 2 (BHP-RATE) - Well-by-Well Validation Errors ---\n');
    fprintf('%-6s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
            'Well', 'BHP_Init(%)', 'BHP_Opt(%)', 'qOs_Init(%)', 'qOs_Opt(%)', 'qWs_Init(%)', 'qWs_Opt(%)');
    fprintf('%-6s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
            '----', '-----------', '-----------', '-----------', '-----------', '-----------', '-----------');
    
    for well_idx = 1:num_wells
        fprintf('%-6d %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f\n', ...
                well_idx, val_bhp_error_old2(well_idx), val_bhp_error_opt2(well_idx), ...
                val_qOs_error_old2(well_idx), val_qOs_error_opt2(well_idx), ...
                val_qWs_error_old2(well_idx), val_qWs_error_opt2(well_idx));
    end
    
else
    %% Non-Nash Case - Show only the chosen model
    
    if playerChoice == 1
        model_name = 'MODEL 1 (RATE-BHP)';
        ref_data = ref_data;
        bhp_error = bhp_error;
        qOs_error = qOs_error;
        qWs_error = qWs_error;
        val_bhp_error_old = val_bhp_error_old;
        val_bhp_error_opt = val_bhp_error_opt;
        val_qOs_error_old = val_qOs_error_old;
        val_qOs_error_opt = val_qOs_error_opt;
        val_qWs_error_old = val_qWs_error_old;
        val_qWs_error_opt = val_qWs_error_opt;
    else
        model_name = 'MODEL 2 (BHP-RATE)';
        ref_data = ref_data;
        bhp_error = bhp_error;
        qOs_error = qOs_error;
        qWs_error = qWs_error;
        val_bhp_error_old = val_bhp_error_old;
        val_bhp_error_opt = val_bhp_error_opt;
        val_qOs_error_old = val_qOs_error_old;
        val_qOs_error_opt = val_qOs_error_opt;
        val_qWs_error_old = val_qWs_error_old;
        val_qWs_error_opt = val_qWs_error_opt;
    end
    
    fprintf('\n--- %s - Well-by-Well Training Errors ---\n', model_name);
    fprintf('%-6s %-10s %-10s %-10s %-15s %-15s %-15s\n', ...
            'Well', 'BHP(%)', 'qOs(%)', 'qWs(%)', 'BHP_ref', 'qOs_ref', 'qWs_ref');
    fprintf('%-6s %-10s %-10s %-10s %-15s %-15s %-15s\n', ...
            '----', '------', '------', '------', '------', '------', '------');
    
    num_wells = size(ref_data.bhp, 2);
    for well_idx = 1:num_wells
        int_bhp = sum(abs(ref_data.bhp(:, well_idx)) .* dt)/barsa;
        int_qOs = sum(abs(ref_data.qOs(:, well_idx)) .* dt);
        int_qWs = sum(abs(ref_data.qWs(:, well_idx)) .* dt);
        
        fprintf('%-6d %-10.2f %-10.2f %-10.2f %-15.1f %-15.2e %-15.2e\n', ...
                well_idx, bhp_error(well_idx), qOs_error(well_idx), qWs_error(well_idx), ...
                int_bhp, int_qOs, int_qWs);
    end
    
    % Validation errors table
    fprintf('\n--- %s - Well-by-Well Validation Errors ---\n', model_name);
    fprintf('%-6s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
            'Well', 'BHP_Init(%)', 'BHP_Opt(%)', 'qOs_Init(%)', 'qOs_Opt(%)', 'qWs_Init(%)', 'qWs_Opt(%)');
    fprintf('%-6s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
            '----', '-----------', '-----------', '-----------', '-----------', '-----------', '-----------');
    
    for well_idx = 1:num_wells
        fprintf('%-6d %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f\n', ...
                well_idx, val_bhp_error_old(well_idx), val_bhp_error_opt(well_idx), ...
                val_qOs_error_old(well_idx), val_qOs_error_opt(well_idx), ...
                val_qWs_error_old(well_idx), val_qWs_error_opt(well_idx));
    end
end

%% Final summary table with model recommendations
fprintf('\n\n=== FINAL MODEL RECOMMENDATIONS ===\n');

if NashAlgorithm
    fprintf('Both models are available for use:\n');
    fprintf('1. MODEL 1 (RATE-BHP): Use for RATE-controlled injection and BHP-controlled production\n');
    fprintf('   - Training Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
            total_bhp_error1, total_qOs_error1, total_qWs_error1);
    fprintf('   - Validation Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
            val_total_bhp_error_opt1, val_total_qOs_error_opt1, val_total_qWs_error_opt1);
    
    fprintf('2. MODEL 2 (BHP-RATE): Use for BHP-controlled injection and RATE-controlled production\n');
    fprintf('   - Training Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
            total_bhp_error2, total_qOs_error2, total_qWs_error2);
    fprintf('   - Validation Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
            val_total_bhp_error_opt2, val_total_qOs_error_opt2, val_total_qWs_error_opt2);
    
    fprintf('\nRECOMMENDATION: Use the appropriate model based on your well control strategy.\n');
    
else
    if playerChoice == 1
        fprintf('Available model: MODEL 1 (RATE-BHP)\n');
        fprintf('Use for RATE-controlled injection and BHP-controlled production\n');
        fprintf('Training Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
                total_bhp_error, total_qOs_error, total_qWs_error);
        fprintf('Validation Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
                val_total_bhp_error_opt, val_total_qOs_error_opt, val_total_qWs_error_opt);
    else
        fprintf('Available model: MODEL 2 (BHP-RATE)\n');
        fprintf('Use for BHP-controlled injection and RATE-controlled production\n');
        fprintf('Training Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
                total_bhp_error, total_qOs_error, total_qWs_error);
        fprintf('Validation Total Errors: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
                val_total_bhp_error_opt, val_total_qOs_error_opt, val_total_qWs_error_opt);
    end
end

fprintf('\n=== CALCULATION METHODOLOGY ===\n');
fprintf('Error = 100 × ∑|optimized - reference| × Δt ÷ max(∑|reference| × Δt, threshold)\n');
fprintf('This gives true time-weighted percentage errors with zero-division protection.\n');
fprintf('BHP threshold: %.1f bar, Rate threshold: %.2e m³/day\n', bhp_threshold/barsa, rate_threshold*(1/day));
fprintf('Zero reference values are handled by setting error to 0%%.\n');



if NashAlgorithm
%% Optimal Model Selection Table - Based on Minimal Validation Errors
fprintf('\n\n=== OPTIMAL MODEL SELECTION BY WELL AND CONTROL TYPE ===\n');
fprintf('Recommendation: Use the model with lowest error for each well/control combination\n\n');

fprintf('%-6s %-12s %-12s %-12s %-20s\n', ...
        'Well', 'Best BHP', 'Best qOs', 'Best qWs', 'Recommended Model');
fprintf('%-6s %-12s %-12s %-12s %-20s\n', ...
        '----', '---------', '---------', '---------', '------------------');

num_wells = 12; % Based on your data
for well_idx = 1:num_wells
    % Determine best model for each variable
    if val_bhp_error_opt1(well_idx) <= val_bhp_error_opt2(well_idx)
        best_bhp = 'Model 1';
        best_bhp_error = val_bhp_error_opt1(well_idx);
    else
        best_bhp = 'Model 2';
        best_bhp_error = val_bhp_error_opt2(well_idx);
    end
    
    if val_qOs_error_opt1(well_idx) <= val_qOs_error_opt2(well_idx)
        best_qOs = 'Model 1';
        best_qOs_error = val_qOs_error_opt1(well_idx);
    else
        best_qOs = 'Model 2';
        best_qOs_error = val_qOs_error_opt2(well_idx);
    end
    
    if val_qWs_error_opt1(well_idx) <= val_qWs_error_opt2(well_idx)
        best_qWs = 'Model 1';
        best_qWs_error = val_qWs_error_opt1(well_idx);
    else
        best_qWs = 'Model 2';
        best_qWs_error = val_qWs_error_opt2(well_idx);
    end
    
    % Determine overall recommended model
    model1_score = (val_bhp_error_opt1(well_idx) <= val_bhp_error_opt2(well_idx)) + ...
                  (val_qOs_error_opt1(well_idx) <= val_qOs_error_opt2(well_idx)) + ...
                  (val_qWs_error_opt1(well_idx) <= val_qWs_error_opt2(well_idx));
    
    model2_score = 3 - model1_score;
    
    if model1_score > model2_score
        recommended_model = 'Model 1 (RATE-BHP)';
    elseif model2_score > model1_score
        recommended_model = 'Model 2 (BHP-RATE)';
    else
        % Tie-breaker: use model with lower average error
        avg_error1 = mean([val_bhp_error_opt1(well_idx), val_qOs_error_opt1(well_idx), val_qWs_error_opt1(well_idx)]);
        avg_error2 = mean([val_bhp_error_opt2(well_idx), val_qOs_error_opt2(well_idx), val_qWs_error_opt2(well_idx)]);
        if avg_error1 <= avg_error2
            recommended_model = 'Model 1 (RATE-BHP)*';
        else
            recommended_model = 'Model 2 (BHP-RATE)*';
        end
    end
    
    fprintf('%-6d %-12s %-12s %-12s %-20s\n', ...
            well_idx, best_bhp, best_qOs, best_qWs, recommended_model);
end

%% Detailed comparison table with actual error values
fprintf('\n--- DETAILED COMPARISON WITH ERROR VALUES ---\n');
fprintf('%-6s %-15s %-15s %-15s %-15s %-15s %-15s\n', ...
    'Well', 'BHP M1(%)', 'BHP M2(%)', 'qOs M1(%)', 'qOs M2(%)', 'qWs M1(%)', 'qWs M2(%)');
fprintf('%-6s %-15s %-15s %-15s %-15s %-15s %-15s\n', ...
    '----', '---------', '---------', '---------', '---------', '---------', '---------');

for well_idx = 1:num_wells
    fprintf('%-6d %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f\n', ...
        well_idx, ...
        val_bhp_error_opt1(well_idx), val_bhp_error_opt2(well_idx), ...
        val_qOs_error_opt1(well_idx), val_qOs_error_opt2(well_idx), ...
        val_qWs_error_opt1(well_idx), val_qWs_error_opt2(well_idx));
end

%% Summary statistics
fprintf('\n--- MODEL PERFORMANCE SUMMARY ---\n');

% Count wins for each model per variable type
bhp_wins_m1 = sum(val_bhp_error_opt1 <= val_bhp_error_opt2);
bhp_wins_m2 = num_wells - bhp_wins_m1;

qOs_wins_m1 = sum(val_qOs_error_opt1 <= val_qOs_error_opt2);
qOs_wins_m2 = num_wells - qOs_wins_m1;

qWs_wins_m1 = sum(val_qWs_error_opt1 <= val_qWs_error_opt2);
qWs_wins_m2 = num_wells - qWs_wins_m1;

fprintf('BHP Performance: Model 1 wins in %d wells, Model 2 wins in %d wells\n', bhp_wins_m1, bhp_wins_m2);
fprintf('qOs Performance: Model 1 wins in %d wells, Model 2 wins in %d wells\n', qOs_wins_m1, qOs_wins_m2);
fprintf('qWs Performance: Model 1 wins in %d wells, Model 2 wins in %d wells\n', qWs_wins_m1, qWs_wins_m2);

% Calculate average errors
avg_bhp_m1 = mean(val_bhp_error_opt1);
avg_bhp_m2 = mean(val_bhp_error_opt2);
avg_qOs_m1 = mean(val_qOs_error_opt1);
avg_qOs_m2 = mean(val_qOs_error_opt2);
avg_qWs_m1 = mean(val_qWs_error_opt1);
avg_qWs_m2 = mean(val_qWs_error_opt2);

fprintf('\nAverage Errors:\n');
fprintf('Model 1 - BHP: %.3f%%, qOs: %.3f%%, qWs: %.3f%%\n', avg_bhp_m1, avg_qOs_m1, avg_qWs_m1);
fprintf('Model 2 - BHP: %.3f%%, qOs: %.3f%%, qWs: %.3f%%\n', avg_bhp_m2, avg_qOs_m2, avg_qWs_m2);

%% Final recommendations by well type
fprintf('\n=== FINAL RECOMMENDATIONS BY WELL TYPE ===\n');

% Identify well types based on patterns
injection_wells = [1:8];  % Wells with zero qOs errors (injection wells)
production_wells = [9:12]; % Wells with non-zero qOs errors (production wells)

fprintf('Injection Wells (1-8):\n');
fprintf('  - BHP control: %s (avg: %.3f%% vs %.3f%%)\n', ...
    iff(avg_bhp_m1 <= avg_bhp_m2 , 'Model 1' , 'Model 2'), avg_bhp_m1, avg_bhp_m2);
fprintf('  - qWs control: %s (avg: %.3f%% vs %.3f%%)\n', ...
    iff(mean(val_qWs_error_opt1(injection_wells)) <= mean(val_qWs_error_opt2(injection_wells)), 'Model 1' , 'Model 2'), ...
    mean(val_qWs_error_opt1(injection_wells)), mean(val_qWs_error_opt2(injection_wells)));

fprintf('Production Wells (9-12):\n');
fprintf('  - BHP control: %s (avg: %.3f%% vs %.3f%%)\n', ...
    iff(mean(val_bhp_error_opt1(production_wells)) <= mean(val_bhp_error_opt2(production_wells)) , 'Model 1' , 'Model 2'), ...
    mean(val_bhp_error_opt1(production_wells)), mean(val_bhp_error_opt2(production_wells)));
fprintf('  - qOs control: %s (avg: %.3f%% vs %.3f%%)\n', ...
    iff(mean(val_qOs_error_opt1(production_wells)) <= mean(val_qOs_error_opt2(production_wells)) , 'Model 1' , 'Model 2'), ...
    mean(val_qOs_error_opt1(production_wells)), mean(val_qOs_error_opt2(production_wells)));
fprintf('  - qWs control: %s (avg: %.3f%% vs %.3f%%)\n', ...
    iff(mean(val_qWs_error_opt1(production_wells)) <= mean(val_qWs_error_opt2(production_wells)) , 'Model 1' , 'Model 2'), ...
    mean(val_qWs_error_opt1(production_wells)), mean(val_qWs_error_opt2(production_wells)));

fprintf('\nOVERALL RECOMMENDATION:\n');
if (bhp_wins_m1 + qOs_wins_m1 + qWs_wins_m1) > (bhp_wins_m2 + qOs_wins_m2 + qWs_wins_m2)
    fprintf('Use Model 1 (RATE-BHP) for most scenarios\n');
elseif (bhp_wins_m1 + qOs_wins_m1 + qWs_wins_m1) < (bhp_wins_m2 + qOs_wins_m2 + qWs_wins_m2)
    fprintf('Use Model 2 (BHP-RATE) for most scenarios\n');
else
    fprintf('Both models perform similarly - choose based on specific well control strategy\n');
end
end
fprintf('\n=== END OF ANALYSIS ===\n');


function result = iff(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end