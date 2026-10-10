
function timesteps = getReportMinisteps(report, schedule)
% Get the timesteps used for the ministeps of a report
% Optional schedule input to verify against control steps
    
    timesteps = [];
    if iscell(report)
        report = struct('ControlstepReports', {report});
    end
    
    % Get all converged ministeps
    for i = 1:numel(report.ControlstepReports)
        cr = report.ControlstepReports{i};
        steprep = cr.StepReports(cellfun(@(x) x.Converged, cr.StepReports));
        current_steps = cellfun(@(x) x.Timestep, steprep);
        timesteps = [timesteps; current_steps]; %#ok
    end
    
    % If schedule is provided, try to match control step count
    if nargin > 1 && ~isempty(schedule)
        expected_steps = numel(schedule.step.val);
        actual_steps = numel(timesteps);
        
        if actual_steps > expected_steps
            % Likely due to timestep cutting - sum ministeps to match control steps
            control_step_times = cumsum(schedule.step.val);
            simulation_times = cumsum(timesteps);
            
            % Find indices where simulation times match control step times
            match_indices = [];
            for t = 1:numel(control_step_times)
                [~, idx] = min(abs(simulation_times - control_step_times(t)));
                match_indices = [match_indices; idx]; %#ok
            end
            
            % Extract timesteps that complete each control step
            if numel(match_indices) == expected_steps
                timesteps = timesteps(match_indices);
            else
                warning('Could not perfectly match control steps. Using all ministeps.');
            end
        end
    end
end