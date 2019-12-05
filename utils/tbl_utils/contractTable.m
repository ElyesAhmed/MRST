function [prodAB, restbl] = contractTable(Acell, Bcell, fds)
%
%
% SYNOPSIS:
%   function [prodAB, restbl] = contractTable(Acell, Bcell, fds)
%
% DESCRIPTION: Compute the contraction of two tensors along the fields given
% by fds
%
% PARAMETERS:
%   Acell - First tensor: Acell{1} gives the tensor in vector form which
%   follows the indexing table given in Acell{2}
%   Bcell - Second tensor: Same structure as first tensor
%   fds   - cell of fields along which the reduction will be made
%
% RETURNS:
%   prodAB - Result of the contraction of the two tensor, as a vector indexed
%   according to indexing table prodAB
%   restbl - Indexing table for prodAB
%
% SEE ALSO:
%   `setupTableMapping`.

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

    A    = Acell{1};
    tblA = Acell{2};
    B    = Bcell{1};
    tblB = Bcell{2};
    
    fds1 = fds{1};
    fds2 = fds{2};
    fdscross = fds{3};
    
    afds1 = fieldnames(tblA);
    afds2 = fieldnames(tblB);
    
    %% sanity checks
    % first check: we should have fds1 (fds2) equal to ofds1 (ofds2).
    ofds1 = afds1;
    fdsToRemove = {fdscross{:}, 'num'};
    for ifield = 1 : numel(fdsToRemove)
        ofds1 = ofds1(~strcmp(ofds1, fdsToRemove{ifield}));
    end
    
    ofds2 = afds2;
    fdsToRemove = {fdscross{:}, 'num'};
    for ifield = 1 : numel(fdsToRemove)
        ofds2 = ofds2(~strcmp(ofds2, fdsToRemove{ifield}));
    end
    
    isdiff = (numel(setdiff(fds1, ofds1)) > 0) | ...
             (numel(setdiff(ofds1, fds1)) > 0) | ...
             (numel(setdiff(fds2, ofds2)) > 0) | ...
             (numel(setdiff(ofds2, fds2)) > 0);
    assert(~isdiff, 'mismatch in table matrix multipication');
    % second check: we do not support for the moment when fds1 and fds2 have
    % common field names (in this case they should be given different names,
    % which we do not do automatically).
    
    assert(numel(intersect(ofds1, ofds2)) == 0, ['repeated names in fields are not ' ...
                        'supported']);
    
    % end of sanity checks
    
    [~, prodmattbl] = setupTableMapping(tblA, tblB, fdscross);
    map1 = setupTableMapping(tblA, prodmattbl, {ofds1{:}, fdscross{:}});
    map2 = setupTableMapping(tblB, prodmattbl, {ofds2{:}, fdscross{:}});
    resfds = {fds1{:}, fds2{:}};
    restbl = projTable(prodmattbl, resfds);
    reducemap = setupTableMapping(prodmattbl, restbl, resfds);
    
    prodAB = reducemap*((map1*A).*(map2*B));
end
