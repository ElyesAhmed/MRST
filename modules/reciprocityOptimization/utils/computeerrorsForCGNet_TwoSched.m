%% ========== Time‑Weighted Error Analysis (adapted for Egg‑model) ==========
% Assumes NashAlgorithm = true and all variables from the Egg script exist.

% Use the common time step vector (same for all schedules)
dt = dt;  % already defined earlier

% Function to extract well data from cell array of well solutions
getWellData = @(ws) struct(...
    'bhp', cell2mat(cellfun(@(x) [x.bhp], ws, 'UniformOutput', false)), ...
    'qOs', cell2mat(cellfun(@(x) [x.qOs], ws, 'UniformOutput', false)), ...
    'qWs', cell2mat(cellfun(@(x) [x.qWs], ws, 'UniformOutput', false)) ...
);

total_time = sum(dt);
bhp_threshold = 1 * barsa;
rate_threshold = 1e-6 * (1/day);

calculateWellError = @(opt_vals, ref_vals, dt, threshold) ...
    100 * sum(abs(opt_vals - ref_vals) .* dt, 1) ./ ...
    max(sum(abs(ref_vals) .* dt, 1), threshold * total_time);

calculateNormWellError = @(opt_vals, ref_vals, dt, threshold) ...
    100 * sum(abs(opt_vals - ref_vals) .* dt, 1) ./ ...
    max(range(ref_vals, 1) * total_time, threshold * total_time);

calculateTotalError = @(opt_vals, ref_vals, dt, threshold) ...
    100 * sum(abs(opt_vals - ref_vals) .* dt, 'all') ./ ...
    max(sum(abs(ref_vals) .* dt, 'all'), threshold * total_time * size(ref_vals, 2));

calculateNormTotalError = @(opt_vals, ref_vals, dt, threshold) ...
    100 * sum(abs(opt_vals - ref_vals) .* dt, 'all') ./ ...
    max(range(ref_vals, 'all') * total_time * size(ref_vals, 2), threshold);

% --- Training errors (Player1 & Player2) ---
ref_data1 = getWellData(wsRef1);
opt_data1 = getWellData(wellSols_opt1);
ref_data2 = getWellData(wsRef2);
opt_data2 = getWellData(wellSols_opt2);

bhp_error1 = calculateNormWellError(opt_data1.bhp, ref_data1.bhp, dt, bhp_threshold);
qOs_error1 = calculateWellError(opt_data1.qOs, ref_data1.qOs, dt, rate_threshold);
qWs_error1 = calculateWellError(opt_data1.qWs, ref_data1.qWs, dt, rate_threshold);

bhp_error2 = calculateNormWellError(opt_data2.bhp, ref_data2.bhp, dt, bhp_threshold);
qOs_error2 = calculateWellError(opt_data2.qOs, ref_data2.qOs, dt, rate_threshold);
qWs_error2 = calculateWellError(opt_data2.qWs, ref_data2.qWs, dt, rate_threshold);

total_bhp_error1 = calculateNormTotalError(opt_data1.bhp, ref_data1.bhp, dt, bhp_threshold);
total_qOs_error1 = calculateTotalError(opt_data1.qOs, ref_data1.qOs, dt, rate_threshold);
total_qWs_error1 = calculateTotalError(opt_data1.qWs, ref_data1.qWs, dt, rate_threshold);

total_bhp_error2 = calculateNormTotalError(opt_data2.bhp, ref_data2.bhp, dt, bhp_threshold);
total_qOs_error2 = calculateTotalError(opt_data2.qOs, ref_data2.qOs, dt, rate_threshold);
total_qWs_error2 = calculateTotalError(opt_data2.qWs, ref_data2.qWs, dt, rate_threshold);

% Protection against Inf/NaN
for f = {'bhp_error1','qOs_error1','qWs_error1','bhp_error2','qOs_error2','qWs_error2'}
    val = eval(f{1});
    val(isinf(val) | isnan(val)) = 0;
    eval([f{1} ' = val;']);
end

fprintf('\n=== TRAINING ERRORS (Time‑Weighted) ===\n');
fprintf('Player 1 (RATE‑BHP) ‑ Total: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
    total_bhp_error1, total_qOs_error1, total_qWs_error1);
fprintf('Player 2 (BHP‑RATE) ‑ Total: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
    total_bhp_error2, total_qOs_error2, total_qWs_error2);

% --- Validation errors (Player1 & Player2) ---
val_ref_data1 = getWellData(ws_valid1);
val_coarse_old_data1 = getWellData(ws_coarse_valid1_old);
val_coarse_opt_data1 = getWellData(ws_coarse_valid1_opt);

val_ref_data2 = getWellData(ws_valid2);
val_coarse_old_data2 = getWellData(ws_coarse_valid2_old);
val_coarse_opt_data2 = getWellData(ws_coarse_valid2_opt);

val_bhp_error_old1 = calculateNormWellError(val_coarse_old_data1.bhp, val_ref_data1.bhp, dt, bhp_threshold);
val_qOs_error_old1 = calculateWellError(val_coarse_old_data1.qOs, val_ref_data1.qOs, dt, rate_threshold);
val_qWs_error_old1 = calculateWellError(val_coarse_old_data1.qWs, val_ref_data1.qWs, dt, rate_threshold);

val_bhp_error_opt1 = calculateNormWellError(val_coarse_opt_data1.bhp, val_ref_data1.bhp, dt, bhp_threshold);
val_qOs_error_opt1 = calculateWellError(val_coarse_opt_data1.qOs, val_ref_data1.qOs, dt, rate_threshold);
val_qWs_error_opt1 = calculateWellError(val_coarse_opt_data1.qWs, val_ref_data1.qWs, dt, rate_threshold);

val_bhp_error_old2 = calculateNormWellError(val_coarse_old_data2.bhp, val_ref_data2.bhp, dt, bhp_threshold);
val_qOs_error_old2 = calculateWellError(val_coarse_old_data2.qOs, val_ref_data2.qOs, dt, rate_threshold);
val_qWs_error_old2 = calculateWellError(val_coarse_old_data2.qWs, val_ref_data2.qWs, dt, rate_threshold);

val_bhp_error_opt2 = calculateNormWellError(val_coarse_opt_data2.bhp, val_ref_data2.bhp, dt, bhp_threshold);
val_qOs_error_opt2 = calculateWellError(val_coarse_opt_data2.qOs, val_ref_data2.qOs, dt, rate_threshold);
val_qWs_error_opt2 = calculateWellError(val_coarse_opt_data2.qWs, val_ref_data2.qWs, dt, rate_threshold);

% Total validation errors
val_total_bhp_error_opt1 = calculateNormTotalError(val_coarse_opt_data1.bhp, val_ref_data1.bhp, dt, bhp_threshold);
val_total_qOs_error_opt1 = calculateTotalError(val_coarse_opt_data1.qOs, val_ref_data1.qOs, dt, rate_threshold);
val_total_qWs_error_opt1 = calculateTotalError(val_coarse_opt_data1.qWs, val_ref_data1.qWs, dt, rate_threshold);

val_total_bhp_error_opt2 = calculateNormTotalError(val_coarse_opt_data2.bhp, val_ref_data2.bhp, dt, bhp_threshold);
val_total_qOs_error_opt2 = calculateTotalError(val_coarse_opt_data2.qOs, val_ref_data2.qOs, dt, rate_threshold);
val_total_qWs_error_opt2 = calculateTotalError(val_coarse_opt_data2.qWs, val_ref_data2.qWs, dt, rate_threshold);

% Protection
for f = {'val_bhp_error_old1','val_qOs_error_old1','val_qWs_error_old1',...
         'val_bhp_error_opt1','val_qOs_error_opt1','val_qWs_error_opt1',...
         'val_bhp_error_old2','val_qOs_error_old2','val_qWs_error_old2',...
         'val_bhp_error_opt2','val_qOs_error_opt2','val_qWs_error_opt2'}
    val = eval(f{1});
    val(isinf(val) | isnan(val)) = 0;
    eval([f{1} ' = val;']);
end

fprintf('\n=== VALIDATION ERRORS ===\n');
fprintf('Player 1 ‑ Initial Coarse : BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
    sum(val_bhp_error_old1)/numel(val_bhp_error_old1), ...
    sum(val_qOs_error_old1)/numel(val_qOs_error_old1), ...
    sum(val_qWs_error_old1)/numel(val_qWs_error_old1));
fprintf('Player 1 ‑ Optimized Coarse: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
    sum(val_bhp_error_opt1)/numel(val_bhp_error_opt1), ...
    sum(val_qOs_error_opt1)/numel(val_qOs_error_opt1), ...
    sum(val_qWs_error_opt1)/numel(val_qWs_error_opt1));
fprintf('Player 2 ‑ Initial Coarse : BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
    sum(val_bhp_error_old2)/numel(val_bhp_error_old2), ...
    sum(val_qOs_error_old2)/numel(val_qOs_error_old2), ...
    sum(val_qWs_error_old2)/numel(val_qWs_error_old2));
fprintf('Player 2 ‑ Optimized Coarse: BHP: %.2f%%, qOs: %.2f%%, qWs: %.2f%%\n', ...
    sum(val_bhp_error_opt2)/numel(val_bhp_error_opt2), ...
    sum(val_qOs_error_opt2)/numel(val_qOs_error_opt2), ...
    sum(val_qWs_error_opt2)/numel(val_qWs_error_opt2));

% Improvement summary
bhp_imp1 = (sum(val_bhp_error_old1)-sum(val_bhp_error_opt1))/max(sum(val_bhp_error_old1),1e-10)*100;
qOs_imp1 = (sum(val_qOs_error_old1)-sum(val_qOs_error_opt1))/max(sum(val_qOs_error_old1),1e-10)*100;
qWs_imp1 = (sum(val_qWs_error_old1)-sum(val_qWs_error_opt1))/max(sum(val_qWs_error_old1),1e-10)*100;
bhp_imp2 = (sum(val_bhp_error_old2)-sum(val_bhp_error_opt2))/max(sum(val_bhp_error_old2),1e-10)*100;
qOs_imp2 = (sum(val_qOs_error_old2)-sum(val_qOs_error_opt2))/max(sum(val_qOs_error_old2),1e-10)*100;
qWs_imp2 = (sum(val_qWs_error_old2)-sum(val_qWs_error_opt2))/max(sum(val_qWs_error_old2),1e-10)*100;

fprintf('\n=== IMPROVEMENT OVER INITIAL ===\n');
fprintf('Player 1 : BHP: %.1f%%, qOs: %.1f%%, qWs: %.1f%%\n', bhp_imp1, qOs_imp1, qWs_imp1);
fprintf('Player 2 : BHP: %.1f%%, qOs: %.1f%%, qWs: %.1f%%\n', bhp_imp2, qOs_imp2, qWs_imp2);

% (Optional) Well‑by‑well table – if needed, adapt using the same variable substitutions.