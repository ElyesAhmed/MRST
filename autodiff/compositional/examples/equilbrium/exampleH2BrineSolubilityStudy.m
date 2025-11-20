%% MRST Example: Solubility Tables for H2-Brine System with RK-EoS vs SW-EoS
%% vs PR-EoS 
%
% This example demonstrates the setup and computation of solubility tables
% for hydrogen (H2) in a brine system using the Redlich-Kwong (RK) Equation
% of State. It also demonstrates the computation of H2 solubility using the
% Peng-Robinson (PR) and the Soreide-Whitson (SW) EoS models with flash calculations.
%
%   Ahmed, E., et al. (2024). Phase behavior and black-oil simulations of
%   hydrogen storage in saline aquifers. Advances in Water Resources, 191, 104772.
%   Ahmed, E., et al. (2025). Modeling and simulation of coupled biochemical 
%   and two-phase compositional flow in underground hydrogen storage
%
% SEE ALSO:
%   `generateComponentTable`, `generateSolubilityTable
%
%
%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

MRST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST. If not, see <http://www.gnu.org/licenses/>.
%}

clear; clc;

mrstModule add biochemistry compositional ad-blackoil ad-core ad-props h2store mrst-gui


% Input Parameters for Temperature, Pressure, and Salinity
min_temp     = 1;                 % [°C]
max_temp     = 99;                % [°C]
min_pressure = 6 * mega * Pascal; % [Pa]
max_pressure = 20 * mega * Pascal;% [Pa]
nbp          = 15;               % Number of pressure points
nbt          = 15;              % Number of temperature points
ms           = 2;                 % Salt molality [mol/kg]
outputDisplay = true;            % Set to true to display generated tables
recompute = true;                 % recompuete solubility  tables

%% Notice on Computational Cost
warning('ComputationalCost:Medium', ...
    ['Please be advised that for large nbp and nbt this example often takes a long time ', ...
    'to run: this script will extract data from https://webbook.nist.gov']);

% Define target output directory and create if it doesn't exist
outputPath = fullfile(mrstOutputDirectory(), 'UHS_PVT', 'H2SolubilityTable');

% This configuration prepares solubility data for the H2-brine system
% under RK-EOS in the context of saline aquifer storage.
% Generate H2 pure Component Table from NIST
pause(0.1);  % Ensure smooth execution between commands
comp_name = 'H2';
disp(['Generating component table for: ', comp_name]);
tab_H2 = generateComponentProperties('min_temp',min_temp, 'max_temp',max_temp, ...
    'n_temp', nbt, 'min_press',min_pressure, 'max_press', max_pressure, ...
    'n_press',nbp, 'comp_name', comp_name,'outputDisplay', outputDisplay, ...
    'outputPath',outputPath);

% We use the Redlich Kwong Eos to obtain the solubility
disp('Generating solubility table for H2-brine mixture...');


tab_sol= generateH2WaterSolubilityTable('min_temp',min_temp, 'max_temp',max_temp, ...
    'n_temp', nbt,'min_press',min_pressure, 'max_press',max_pressure, ...
    'n_press', nbp, 'ms', ms,'outputDisplay', outputDisplay,'outputPath',outputPath, ...
    'reCompute', recompute);


% Load ePC-SAFT Data
epcsaft = load(fullfile(ROOTDIR,'..','modules','h2store','examples','data',...
    'PcSaftSolubilityTable','ePcSaftH2BrineData.mat'));
state = epcsaft.state;

% Define constants
n = nbt;                         % Number of temperature points
indexH2 = 2;                      % Index of H2 component in liquid phase
indexH2O = 1;                   % Index of water in vapor phase

pressure = reshape(tab_sol.("pressure [Pa]"), [], n);
temperature = reshape(tab_sol.("# temperature [°C]"), [], n) + 273.15;  % Convert to Kelvin

%% calculate H2 solubility from PR and SW EoS models (Flash calculations)
%% Define Composition Mixture
% Create a compositional mixture with water and hydrogen components
compFluid = TableCompositionalMixture({'Water', 'Hydrogen'}, {'H2O', 'H2'});
disp(compFluid);

%% Initialize Thermodynamic Model
% Choose an equation of state model for the calculations
eosNamesw = 'sw'; % Soreide-Whitson (SW) model
eosNamepr = 'pr'; % Peng-Robinson (PR) model
eosModelsw = EquationOfStateModel([], compFluid, eosNamesw);
eosModelpr = EquationOfStateModel([], compFluid, eosNamepr);
eosModelsw.msalt=ms;
z0 = [0.8, 0.2]; % Initial composition
pressSWPR=tab_sol.("pressure [Pa]");
tempSWPR=tab_sol.("# temperature [°C]")+ 273.15;

nc = numel(pressSWPR);
namecp = eosModelsw.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indH2O= find(strcmp(namecp,'H2O'));

[Lsw, xsw, ~] = standaloneFlash(pressSWPR, tempSWPR, z0, eosModelsw);
[Lpr, xpr, ~] = standaloneFlash(pressSWPR, tempSWPR, z0, eosModelpr);
xliqH2sw=xsw(:,indH2);
xliqH2pr=xpr(:,indH2);


%% Solubility of H2 in Brine as a Function of Pressure and Temperature
figure;
hold on;
plotColors = {'k', 'r', 'b', 'm'}; % Colors for each pressure level
kk = 1:int16((n - 1) / 4):n;
legends = cell(16,1); % Initialize legends array
% Loop over the four selected pressure levels
for j = 1:4
    i = kk(j);
    % Convert the pressure value to MPa and set up legend label
    P_val = unique(pressSWPR(i:n:end));
    P_legend = [num2str(int16(convertTo(P_val, mega * Pascal)))];
    % Plot SW results with solid lines
    plot(tempSWPR(i:n:n * n), xliqH2sw(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '-', 'LineWidth', 1.5);
    legends{4*j-3} = ['SW EoS ' P_legend ' MPa']; % Add legend entry
    % Plot RK-EoS results with dashed lines for the same pressures
    plot(tempSWPR(i:n:n * n), tab_sol.x_H2(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '--', 'LineWidth', 1.5);
    legends{4*j-2} = ['RK EoS' P_legend ' MPa']; % Add legend entry
  % Plot PR-EoS results with dot-dash lines for the same pressures
    plot(tempSWPR(i:n:n * n), xliqH2pr(i:n:n * n), ...
       'Color', plotColors{j}, 'LineStyle', '-.', 'LineWidth', 1.5);
    legends{4*j-1} = ['PR EoS' P_legend ' MPa']; % Add legend entry

     plot(state.T(i:n:n * n), state.X_L(i:n:n * n, indexH2), ...
        'Color', plotColors{j}, 'LineStyle', ':', 'LineWidth', 1.5);
    legends{4*j} = ['ePC-SAFT ' P_legend ' MPa']; % Add legend entry   
end

xlabel('Temperature (K)', 'FontSize', 14);
ylabel('x_{H2}', 'FontSize', 14);
legend(legends, 'FontSize', 10, 'Location', 'best');
title('H2 Solubility in salt water: SW-EoS vs. PR-EoS vs. RK-EoS vs ePC-SAFT', 'FontSize', 12);
hold off;







%% Plot the Amount of H2 in Water
figure;
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(tab_sol.x_H2, [], n), 'LineStyle', 'none');
colorbar; % Add a color bar for reference
title('H2 Solubility in salt Water (RK eos)', 'FontSize', 14);
xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Pressure (MPa)', 'FontSize', 14);
set(gca, 'FontSize', 12); % Set font size for axes
grid on; % Add grid for better readability

figure;
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(xliqH2sw, [], n), 'LineStyle', 'none');
colorbar; % Add a color bar for reference
title('H2 Solubility in salt Water (SW eos)', 'FontSize', 14);
xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Pressure (MPa)', 'FontSize', 14);
set(gca, 'FontSize', 12); % Set font size for axes
grid on; % Add grid for better readability

figure;
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(xliqH2pr, [], n), 'LineStyle', 'none');
colorbar; % Add a color bar for reference
title('H2 Solubility in salt Water (PR eos)', 'FontSize', 14);
xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Pressure (MPa)', 'FontSize', 14);
set(gca, 'FontSize', 12); % Set font size for axes
grid on; % Add grid for better readability



% Force recompute tables
if ~recompute
    disp('You have changed P, T, or ms! Make sure "recompute" is set to true to recalculate tables.');
end
