function [data_table, status, file_path] = generatePureComponentTable(min_temp, max_temp, nbt, comp_name, output_dir, varargin)
    % generateComponentTable Generates a pure component table for the specified component.
    %
    % This function creates a thermodynamic table by running an external Python
    % script with specified parameters for temperature, pressure, and sampling points.
    % The generated file is saved in the specified directory, and the table is
    % read into MATLAB as a table output.
    %
    % Inputs:
    %   min_temp       - Minimum temperature in Celsius.
    %   max_temp       - Maximum temperature in Celsius.
    %   nbt            - Number of temperature sampling points.
    %   comp_name      - Name of the component (e.g., 'H2O').
    %   output_dir     - Directory where the generated file will be saved.
    %
    % Outputs:
    %   data_table     - MATLAB table containing the generated data.
    %   status         - Status of the system command execution (0 for success).
    %   file_path      - Full path of the generated CSV file.

    % SEE ALSO:
    %   generateSolubilityTable, generateComponentTable.

    %{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opts.outputDisplay = true;  % Default to display output
    opts = merge_options(opts,varargin{:});
    if ~isfolder(output_dir)
        mkdir(output_dir); % Create the directory if it doesn't exist

    end

    % Construct the filename based on the defined parameters
    file_name = sprintf('%spurevalues_%.1f_to_%.1f_C.csv', lower(comp_name), ...
        min_temp, max_temp);

    % Full file path
    file_path = fullfile(output_dir, file_name);

    % Define the command to run the Python script, passing parameters to control output
    command = sprintf('python3 /python_scripts/make_purecomponent_table_H2.py -t1 %.1f -t2 %.1f -nt %d -c %s -o %s', ...
        min_temp, max_temp, nbt, lower(comp_name), file_path);

    % Print the command for debugging
    disp(['Executing command: ', command]);

    % Execute the command using the system function
    [status, result] = system(command);

    % Check the status and display the result
    if status == 0
        disp('Command executed successfully:');
        disp(result);
    else
        disp('Error occurred while executing the command:');
        disp(result);
        return; % Exit the function if the command fails
    end

    % Display path information for verification
    fprintf('File will be saved to: %s\n', file_path);

    % Read the CSV file into a table
    try
        data_table = readtable(file_path);
        disp('Table read successfully:');
        if opts.outputDisplay
            disp(data_table);
        end
    catch ME
        disp('Error occurred while reading the table:');
        disp(ME.message);
        data_table = []; % Return empty if reading fails
    end
end
