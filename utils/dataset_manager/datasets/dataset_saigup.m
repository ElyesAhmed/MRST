function [info, present] = dataset_saigup()
% Info function for SAIGUP dataset. Use getDatasetInfo or getAvailableDatasets for practical purposes.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
    [info, present] = datasetInfoStruct(...
        'name', 'SAIGUP', ...
        'website', 'https://www.sintef.no/projectweb/mrst/downloadable-resources/download-saigup-data-set/', ...
        'fileurl', 'https://www.sintef.no/globalassets/project/mrst/saigup.tar.gz', ...
        'hasGrid', true, ...
        'hasRock', true, ...
        'hasFluid', false, ...
        'cells', 78720, ...
        'examples', {'showSAIGUP', ...
                     'book:showGridSAIGUP', ...
                     'book:showRockSAIGUP', ...
                     'book:interactiveSAIGUP', ...
                     'book:coarsenSAIGUP', ...
                     'coarsegrid:coarsegridExampleSAIGUP', ...
                     'incomp:incompExampleSAIGUP1ph', ...
                     'incomp:incompExampleSAIGUP2ph', ...
                     'diagnostics:diagnostViewerDemo'}, ...
        'description','The SAIGUP project generated a large number of realistic geological models representing shallow marine environments. This model is one specific realization that has faults, inactive cells, disconnected components, but no pinchouts.',... 
        'filesize',    2.535, ...
        'modelType', 'Corner-point' ...
         );
end
