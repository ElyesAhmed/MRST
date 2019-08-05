function [ws, states, reports] = getPackedSimulatorOutput(problem, varargin)
%Get output from a packed simulation problem
%
% SYNOPSIS:
%   [ws, states, reports] = getPackedSimulatorOutput(problem)
%
% REQUIRED PARAMETERS:
%   problem - Problem generated using 'packSimulationProblem'.
%
% OPTIONAL PARAMETERS:
%   readFromDisk - Indicating if states and reporst will be read from disk,
%                  or returned as OutputHandler instances. Reading from
%                  disk can take some time and is recommended when further
%                  changes to output is desired.
%   readWellSolsFromDisk - See above. Applies for wellSols only.
%
% RETURNS:
%   ws      - Well output.
%   states  - States for each simulated step.
%   reports - Report for each of the time-steps.
% 
% EXAMPLE:
%   demoPackedProblems
%
% SEE ALSO:
%   packSimulationProblem, getMultiplePackedSimulatorOutputs

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('readFromDisk', true, 'readWellSolsFromDisk', true);
    opt = merge_options(opt, varargin{:});
    
    nstep = numel(problem.SimulatorSetup.schedule.step.val);
    
    sh = problem.OutputHandlers.states;
    wh = problem.OutputHandlers.wellSols;
    rh = problem.OutputHandlers.reports;
    
    ndata = sh.numelData();
    [ws, states, reports] = deal(cell(ndata, 1));
    wantWells = false;
    for i = 1:numel(problem.SimulatorSetup.schedule.control)
        ctrl = problem.SimulatorSetup.schedule.control(i);
        if isfield(ctrl, 'W') && ~isempty(ctrl.W)
            wantWells = true;
            break
        end
    end

    sn = sprintf('%s (%s)', problem.BaseName, problem.Name);
    if nstep == ndata
        fprintf('Found complete data for %s: %d steps present\n', sn, ndata);
    elseif ndata > nstep
        warning('Found too much data for %s: %d of %d steps present. Case may have been redefined!\n', sn, ndata, nstep);
    elseif ndata > 0
        fprintf('Found partial data for %s: %d of %d steps present\n', sn, ndata, nstep);
    else
        fprintf('Did not find data for %s\n', sn);
    end
    wellOutputMissing = wantWells && wh.numelData() == 0;
    for i = 1:ndata
        if nargout > 1 && opt.readFromDisk
            states{i} = sh{i};
        end
        if wantWells && opt.readWellSolsFromDisk
            if wellOutputMissing
                if isempty(states{i})
                    ws{i} = sh{i}.wellSol;
                else
                    ws{i} = states{i}.wellSol;
                end
            else
                try
                    ws{i} = wh{i};
                catch
                    ws{i} = states{i}.wellSol;
                end
            end
        end
        if nargout > 2 && opt.readFromDisk
            try
                reports{i} = rh{i};
            catch
                reports{i} = [];
            end
        end
    end

    if ~opt.readFromDisk
        % Just return handlers instead
        states = sh;
        reports = rh;
    end
    
    if ~opt.readWellSolsFromDisk
        ws = wh;
    end
end