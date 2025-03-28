function [info, present] = dataset_spe10_layer_compositional()
% Info function for spe10 layer 1 compositional dataset. Use getDatasetInfo or getAvailableDatasets for practical purposes.

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
    [info, present] = datasetInfoStruct(...
        'name'       , 'SPE10_Layer_Compositional', ...
        'website'    , '', ...
        'fileurl'    , 'https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/SPE10_Layer_Compositional.zip', ...
        'hasGrid'    , false, ...
        'hasRock'    , false, ...
        'description', ['Water-alternating gas injection in layer 1 of SPE10 Model 2 ', ...
                        'using a six-component model. Example from Moncorgé et al, '  , ...
                        'J. Comput. Phys, 2018, doi: 10.1016/j.jcp.2018.05.048. '     , ...
                        'See `spe10_layer_compositional`.'                            ], ...
        'hasFluid'   , true, ...
        'examples'   , {'spe10_layer_compositional.m'}, ...
        'filesize'   , 7.6 ...
        );
end
