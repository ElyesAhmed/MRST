 aa =2;
 states = scenarios{aa}.states;
 ws = scenarios{aa}.ws;
 model = scenarios{aa}.model;
% %schedule = scenarios{aa}.schedule;
% --- Compute total H₂ injected (excluding discharge and shut periods) ---
% Define the names of wells/periods to EXCLUDE
exclude_names = {'discharge', 'shut'};  % adjust as needed

% total_injected_H2 = 0;
% for t = 1:length(schedule.step.val)
%     % Get control index for this timestep
%     % ctrl = schedule.step.control(t);
%     % % Get all wells active in this control
%     % wells = schedule.control(ctrl).W;
%     % 
%     % % --- Skip this timestep if any active well has an excluded name ---
%     % skip = false;
%     % for w = 1:numel(wells)
%     %     if any(strcmp(wells(w).name, exclude_names))
%     %         skip = true;
%     %         break;
%     %     end
%     % end
%     % if skip
%     %     continue;
%     % end
%     % 
%     % % --- Now we are in an injection period (no excluded wells) ---
%     % % Find the well(s) that are actually injecting (sign=+1)
%     % % inj_wells = wells([wells.sign] == 1);
%     % % if isempty(inj_wells)
%     % %     continue;  % no injecting well (should not happen)
%     % % end
%     % 
%     % % Get the well solution for this timestep
%     ws_t = ws{t};
% 
%     % For each injecting well, extract H₂ mass rate and accumulate
%     % for iw = 1:numel(inj_wells)
%     %     wname = inj_wells(iw).name;
%     %     % Find this well in the well solution by name
%     %     ws_idx = find(strcmp({ws_t.name}, wname), 1);
%     %     if isempty(ws_idx)
%     %         error('Well "%s" not found in well solution at timestep %d', wname, t);
%     %     end
%     H2_mass_rate = sum(ws_t.ComponentTotalFlux(:,idx_H2));
%     dt = schedule.step.val(t);  % seconds
%     total_injected_H2 = [total_injected_H2; (H2_mass_rate / mc_H2) * dt];
%     %end
% end
% fprintf('Total H₂ injected (excluding discharge and shut periods): %.2f mol\n', total_injected_H2);
eosNames = model.EOSModel.CompositionalMixture.names;
idx_H2 = find(strcmp(eosNames, 'H2'));
mc_H2 = model.EOSModel.CompositionalMixture.molarMass(idx_H2);
% Option B: from well solutions (more accurate)
% NOTE: ComponentTotalFlux is a MASS rate [kg/s] (EquationOfStateComponent
% .getComponentDensity weights the mass phase density by mass fraction),
% so it must be divided by the component molar mass to get mol/s.
total_injected_H2 = 0;
for t = 1:length(schedule.step.val)
    ws_t = ws{t};
    ic = schedule.step.control(t).W;
    idx = strcmp({iW.name}, 'Discharge');
    W_ic = schedule.control(idx).W;
    if ~any(idx)
        H2_mass_rate = ws_t(idx).ComponentTotalFlux(idx_H2);
        dt = schedule.step.val(t);
        total_injected_H2 = total_injected_H2 + (H2_mass_rate / mc_H2) * dt;
    end
end
fprintf('Total injected H₂: %.2f mol\n', total_injected_H2);
% Compute per-cell consumption for each reaction
[H2_rate_meth, H2_cum_meth] = computeH2Consumption(states, schedule, model, 1);
[H2_rate_aceto, H2_cum_aceto] = computeH2Consumption(states, schedule, model, 2);
[H2_rate_srb, H2_cum_srb] = computeH2Consumption(states, schedule, model, 3);
% Sum over all cells to get total cumulative at each timestep
total_cum_meth = sum(H2_cum_meth, 1)';   % column vector
total_cum_aceto = sum(H2_cum_aceto, 1)';
total_cum_srb = sum(H2_cum_srb, 1)';
total_cum_all = total_cum_meth + total_cum_aceto + total_cum_srb;
% Time vector (days)
time_days = cumsum(schedule.step.val) / (24 * 3600);
% Plot
figure;
plot(time_days, total_cum_meth, 'b-', 'LineWidth', 2); hold on;
plot(time_days, total_cum_aceto, 'r-', 'LineWidth', 2);
plot(time_days, total_cum_srb, 'g-', 'LineWidth', 2);
plot(time_days, total_cum_all, 'k--', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Cumulative H₂ consumed (mol)');
legend('Methanogens', 'Acetogens', 'SRB', 'Total', 'Location', 'best');
title('Total H₂ Consumption over Time');
grid on;
% Final cumulative consumption per cell for each reaction
final_cum_meth = H2_cum_meth(:, end);
final_cum_aceto = H2_cum_aceto(:, end);
final_cum_srb = H2_cum_srb(:, end);
final_cum_all = final_cum_meth + final_cum_aceto + final_cum_srb;
% Dimensionless length (cell centers normalized)
x_coords = model.G.cells.centroids(:, 1);
x_norm = (x_coords - min(x_coords)) / (max(x_coords) - min(x_coords));
% Plot each reaction separately
figure;
plot(x_norm, final_cum_meth, 'b-', 'LineWidth', 2); hold on;
plot(x_norm, final_cum_aceto, 'r-', 'LineWidth', 2);
plot(x_norm, final_cum_srb, 'g-', 'LineWidth', 2);
plot(x_norm, final_cum_all, 'k--', 'LineWidth', 2);
xlabel('Dimensionless length (distance from injector)');
ylabel('Cumulative H₂ consumed per cell (mol)');
legend('Methanogens', 'Acetogens', 'SRB', 'Total', 'Location', 'best');
title('Spatial Distribution of H₂ Consumption');
grid on;
sum(final_cum_all)./total_injected_H2*100