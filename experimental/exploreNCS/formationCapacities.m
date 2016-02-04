%% Upper bound trapping capacities: Norwegian Continental Shelf formations

%% Specify the formations to study
% We list the formations explicity since some formations are sealing
% or its grid may not be generated by getAtlasGrid()

names = [getBarentsSeaNames() getNorwegianSeaNames() getNorthSeaNames()];

% Remove certain formation names:
names = names(~strcmpi(names,'Nordmelafm'));
names = names(~strcmpi(names,'Rorfm'));
names = names(~strcmpi(names,'Notfm'));
names = names(~strcmpi(names,'Knurrfm'));       % @@ can be constructed??
names = names(~strcmpi(names,'Fruholmenfm'));   % @@
names = names(~strcmpi(names,'Cookfm'));
names = names(~strcmpi(names,'Dunlingp'));
names = names(~strcmpi(names,'Paleocene'));

%N = 5;
% Load res containing formation names and their coarsening levels.
% Or get res from testing_coarsening_levels()
load coarsening_levels_dx3000meter.mat
shared_names = intersect(names, {res{:,1}});
assert( numel(shared_names) >= numel(names) )
assert( all(strcmpi(sort(shared_names),sort(names)))==1 )

fmCapacities = cell(numel(names),1);
savePlot = true;
figDirName = 'mapPlots_dx3000meters';
mkdir(figDirName);
rhoCref = 760 * kilogram/meter^3;

keep_consistent_scale = true;
width  = 351000; % Bryne's width
height = 420500; % Utsira's height


for i = 1:numel(names)
   
    inx          = find(strcmp(names{i},{res{:,1}}));
    N            = res{inx,2};
    [Gt, rock2D] = getFormationTopGrid(names{i}, N);
    if any(isnan(rock2D.perm))
        rock2D.perm = 500*milli*darcy * ones(Gt.cells.num,1);
        fprintf('\n\nUsing default permeability:\n')
    end
    if any(isnan(rock2D.poro))
        rock2D.poro = 0.25 * ones(Gt.cells.num,1);
        fprintf('Using default porosity:')
    end
    seainfo         = getSeaInfo(names{i}, rhoCref);
    fmCapacities{i} = getTrappingInfo(Gt, rock2D, seainfo, 'mapPlotOn',true, 'fmName',names{i});
    title(names{i})
    
    if keep_consistent_scale
        set(gcf,'Position',[1 1 582 697])
        % resize axis limits to largest formation size:
        cx = min(Gt.cells.centroids(:,1)) + (max(Gt.cells.centroids(:,1)) - min(Gt.cells.centroids(:,1)))/2;
        cy = min(Gt.cells.centroids(:,2)) + (max(Gt.cells.centroids(:,2)) - min(Gt.cells.centroids(:,2)))/2;
        xlim([cx-width/2,  cx+width/2])
        ylim([cy-height/2, cy+height/2])
    end
    
    if savePlot
        title(names{i})
        %title('')
        drawnow
        pause(1)
        if keep_consistent_scale
            export_fig(gcf,[figDirName '/' names{i} '_ref',num2str(N)], '-png','-transparent','-nocrop')
        else
            export_fig(gcf,[figDirName '/' names{i} '_ref',num2str(N) '_scaled'], '-png','-transparent')
        end
        close
    end
end


%% Create table of trapping output

fprintf('\n Name & Structural & Residual & Dissolved & Total in Gt \\\\ \n')
for i = 1:numel(names)
    
    x = fmCapacities{i}.breakdown;
    fprintf('%15s & %8.2f (%5.1f) & %8.2f (%5.1f) & %8.2f (%5.1f) & %8.2f \\\\ \n', names{i}, ...
        x.structural_trapping_capacity, (x.structural_trapping_capacity/x.total_trapping_capacity)*100, ...
        x.residual_trapping_capacity, (x.residual_trapping_capacity/x.total_trapping_capacity)*100, ...
        x.dissolved_trapping_capacity, (x.dissolved_trapping_capacity/x.total_trapping_capacity)*100, ...
        x.total_trapping_capacity );
    
end
    
%% Create table of parameter values used to compute trapping capacities

info = getSeaInfo('NorwegianSea', 760 * kilogram/meter^3);

fprintf('\n Norwegian Sea Parameter Values \n')
fprintf('Parameter                  & Value     & Unit & Reference \\\\ \n')
fprintf('Sea depth                  & %d        & meter & \\\\ \n', getfield(info,'seafloor_depth')*meter)
fprintf('Injection depth            & %d        &       & \\\\ \n', [])
fprintf('Thermal gradient           & %4.1f     & C/km  & \\\\ \n', getfield(info,'temp_gradient'))
fprintf('Seabed temperature         & %d        & C     & \\\\ \n', getfield(info,'seafloor_temp'))
fprintf('Residual water saturation  & %4.2f     &       & \\\\ \n', getfield(info,'res_sat_wat'))
fprintf('Residual {\\co} saturation  & %4.2f     &       & \\\\ \n', getfield(info,'res_sat_co2'))
fprintf('Rock porosity              & %d        &       & \\\\ \n', [])
fprintf('Rock net-to-gross          & %d        &       & \\\\ \n', [])
fprintf('Water density              & %d        & kg/m^3 & \\\\ \n', getfield(info,'water_density')*kilogram/meter^3)
fprintf('{\\co} solubility in brine  & %d        & kg/m^3 & \\\\ \n', getfield(info,'co2_solubility')*kilogram/meter^3)


info = getSeaInfo('BarentsSea', 760 * kilogram/meter^3);

fprintf('\n Barents Sea Parameter Values \n')
fprintf('Parameter                  & Value     & Unit & Reference \\\\ \n')
fprintf('Sea depth                  & %d        & meter & \\\\ \n', getfield(info,'seafloor_depth')*meter)
fprintf('Injection depth            & %d        &       & \\\\ \n', [])
fprintf('Thermal gradient           & %4.1f     & C/km  & \\\\ \n', getfield(info,'temp_gradient'))
fprintf('Seabed temperature         & %d        & C     & \\\\ \n', getfield(info,'seafloor_temp'))
fprintf('Residual water saturation  & %4.2f     &       & \\\\ \n', getfield(info,'res_sat_wat'))
fprintf('Residual {\\co} saturation  & %4.2f     &       & \\\\ \n', getfield(info,'res_sat_co2'))
fprintf('Rock porosity              & %d        &       & \\\\ \n', [])
fprintf('Rock net-to-gross          & %d        &       & \\\\ \n', [])
fprintf('Water density              & %d        & kg/m^3 & \\\\ \n', getfield(info,'water_density')*kilogram/meter^3)
fprintf('{\\co} solubility in brine  & %d        & kg/m^3 & \\\\ \n', getfield(info,'co2_solubility')*kilogram/meter^3)


info = getSeaInfo('NorthSea', 760 * kilogram/meter^3);

fprintf('\n North Sea Parameter Values \n')
fprintf('Parameter                  & Value     & Unit & Reference \\\\ \n')
fprintf('Sea depth                  & %d        & meter & \\\\ \n', getfield(info,'seafloor_depth')*meter)
fprintf('Injection depth            & %d        &       & \\\\ \n', [])
fprintf('Thermal gradient           & %4.1f     & C/km  & \\\\ \n', getfield(info,'temp_gradient'))
fprintf('Seabed temperature         & %d        & C     & \\\\ \n', getfield(info,'seafloor_temp'))
fprintf('Residual water saturation  & %4.2f     &       & \\\\ \n', getfield(info,'res_sat_wat'))
fprintf('Residual {\\co} saturation  & %4.2f     &       & \\\\ \n', getfield(info,'res_sat_co2'))
fprintf('Rock porosity              & %d        &       & \\\\ \n', [])
fprintf('Rock net-to-gross          & %d        &       & \\\\ \n', [])
fprintf('Water density              & %d        & kg/m^3 & \\\\ \n', getfield(info,'water_density')*kilogram/meter^3)
fprintf('{\\co} solubility in brine  & %d        & kg/m^3 & \\\\ \n', getfield(info,'co2_solubility')*kilogram/meter^3)



%% Average rock properties of formations
fprintf('\n Name & Avg Poro (max,min) & Avg Perm (max,min) [mD] & Avg NTG (max,min) \n')
for i = 1:numel(names)
   
    r       = fmCapacities{i}.rock2D;
    r.perm  = convertTo(r.perm, milli*darcy);
    if ~isfield(r,'ntg')
        r.ntg = [];
    end
    fprintf('%16s  & %4.2f (%4.2f, %4.2f)  & %3.0f (%3.0f, %3.0f)  & %4.2f (%4.2f, %4.2f) \\\\ \n', names{i}, ...
        mean(r.poro), max(r.poro), min(r.poro), ...
        mean(r.perm), max(r.perm), min(r.perm), ...
        mean(r.ntg),  max(r.ntg),  min(r.ntg) );
    
end


%% Surface area of formations
fprintf('\n Name & Top surface area [km2] \n')
for i = 1:numel(names)
   
    area = sum(fmCapacities{i}.Gt.cells.volumes);
    area = convertTo(area,(kilo*meter)^2);
    fprintf('%16s  & %8.2f \\\\ \n', names{i}, area);
    
end

%% Number of spill-regions per formation
fprintf('\n Name & Number of Spill-regions \n')
for i = 1:numel(names)

    sptr = fmCapacities{i}.ta.trap_regions;
    num_sptr = max( unique(sptr( sptr > 0 )) );
    fprintf('%16s  & %8.2f \\\\ \n', names{i}, num_sptr);
end

%% Plot Structural Capacity versus coarsening level:
figure(100);
hold on
for i = 1:numel(names)
    
   strap_cap = fmCapacities{i}.breakdown.structural_trapping_capacity;
   c_level = N;
   plot(c_level, strap_cap, 'o')
    
end
legend





