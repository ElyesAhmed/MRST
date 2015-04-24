function [info, present] = dataset_bedmodels1()
% Info function for bedModels1 dataset. Use getDatasetInfo or getAvailableDatasets for practical purposes.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
        'name', 'BedModels1', ...
        'website', 'http://brage.bibsys.no/xmlui/handle/11250/259148', ...
        'fileurl', 'http://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/models_laods.zip', ...
        'description', 'Benchmark models for local upscaling created by Lars Vingli Odsaeter in his thesis "Numerical Aspects of Flow Based Local Upscaling" at the Norwegian University of Science and Technology (2013).', ...
        'hasGrid', true, ...
        'hasRock', true, ...
        'hasFluid', false, ...
        'cells',    21778, ...
        'source', 'Lars Vingli Odsaeter / NTNU', ...
        'filesize',    3.3, ...
        'modelType', 'Corner point' ...
         );
end
