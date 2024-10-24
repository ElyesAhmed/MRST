function [W, trajectory] = addWellFromTrajectory(W, G, rock, traj, varargin)
% This function adds a well based on a piecewise linear trajectory traj.
% All function arguments are the same as for addWell except for in place of
% cellInx, there should be a nx3 matrix of trajectory coordinates.

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

opt = struct('faces',                  [], ...
             'exteriorFaceCorrection', false, ...
             'refDepth',               [], ...
             'errorOnEmptyCellSet',    true);
[opt, other] = merge_options(opt, varargin{:});
assert(G.griddim == 3, 'Function only compatible with 3D grids.');
assert(size(traj,2)==3, 'Trajectory is expected to be nx3 array.');

if ~isfield(G.faces, 'bbox')
    G = addBoundingBoxFields(G);
end

trajectory = computeTraversedCells(G, traj, 'faces', opt.faces, ...
            'exteriorFaceCorrection', opt.exteriorFaceCorrection);

if isempty(trajectory.cell)
    if opt.errorOnEmptyCellSet
        error('Did not find any traversed cells for trajectory.');
    else
        % return shut dummy well
        W = addWell(W, G, rock, 1, other{:}, 'refDepth', opt.refDepth);
        W(end).status  = false;
        W(end).cstatus = false;
    end
else
    if isempty(opt.refDepth) && isfield(G.cells, 'centroids')
        opt.refDepth = G.cells.centroids(trajectory.cell(1), 3);
    end
    % multiply segments by weight (non-unit for segments shared by multiple cells)
    seg = bsxfun(@times, trajectory.vec, trajectory.weight);
    W = addWell(W, G, rock, trajectory.cell, 'lineSegments', seg, 'refDepth', opt.refDepth, other{:});
end
end

