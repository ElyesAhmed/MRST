function [reducedMaps, reducedU] = reduceControlMappings(maps, u, varargin)
% REDUCECONTROLMAPPINGS Reduce/filter control mappings and control vector
%
% SYNOPSIS:
%   [reducedMaps, reducedU] = reduceControlMappings(maps, u)
%   [reducedMaps, reducedU] = reduceControlMappings(maps, u, 'type', {'bhp', 'rate'})
%   [reducedMaps, reducedU] = reduceControlMappings(maps, u, 'wellNo', [1, 3])
%   [reducedMaps, reducedU] = reduceControlMappings(maps, u, 'stepNo', 1)
%   [reducedMaps, reducedU] = reduceControlMappings(maps, u, 'isTarget', true)
%   [reducedMaps, reducedU] = reduceControlMappings(maps, u, 'excludeLims', true)
%   [reducedMaps, reducedU] = reduceControlMappings(maps, u, 'type', 'bhp', 'wellNo', 1, 'excludeLims', true)
%
% PARAMETERS:
%   maps      - Control mapping structure
%   u         - Control vector
%   varargin  - Filter criteria as name-value pairs:
%               'type'        - Cell array of control types to include
%               'wellNo'      - Array of well numbers to include
%               'stepNo'      - Array of step numbers to include
%               'isTarget'    - Logical value to filter by target status
%               'excludeLims' - Logical value to exclude limit controls (default: false)
%
% RETURNS:
%   reducedMaps - Filtered mapping structure
%   reducedU    - Filtered control vector

    if isempty(varargin)
        % No filters specified, return original
        reducedMaps = maps;
        reducedU = u;
        return;
    end
    
    % Parse input arguments
    opt = struct('type', [], 'wellNo', [], 'stepNo', [], 'isTarget', [], 'excludeLims', false);
    opt = merge_options(opt, varargin{:});
    
    % Create logical index for filtering
    nControls = numel(maps.wellNo);
    keepIdx = true(nControls, 1);
    
    % Apply filters
    if ~isempty(opt.type)
        if ischar(opt.type)
            opt.type = {opt.type};
        end
        typeMatch = false(nControls, 1);
        for i = 1:nControls
            typeMatch(i) = any(strcmp(maps.type{i}, opt.type));
        end
        keepIdx = keepIdx & typeMatch;
    end
    
    if ~isempty(opt.wellNo)
        keepIdx = keepIdx & ismember(maps.wellNo, opt.wellNo);
    end
    
    if ~isempty(opt.stepNo)
        keepIdx = keepIdx & ismember(maps.stepNo, opt.stepNo);
    end
    
    if ~isempty(opt.isTarget)
        keepIdx = keepIdx & (maps.isTarget == opt.isTarget);
    end
    
    % Apply excludeLims filter
    if opt.excludeLims
        % Identify which controls are limits (not targets)
        isLimitControl = ~maps.isTarget;
        keepIdx = keepIdx & ~isLimitControl;
    end
    
    % Apply filtering
    reducedU = u(keepIdx);
    
    % Filter maps structure
    reducedMaps = struct();
    fields = fieldnames(maps);
    for i = 1:numel(fields)
        fieldName = fields{i};
        fieldData = maps.(fieldName);
        
        if iscell(fieldData)
            reducedMaps.(fieldName) = fieldData(keepIdx);
        else
            reducedMaps.(fieldName) = fieldData(keepIdx, :);
        end
    end
end

% Helper function for option merging (commonly used in MRST)
function opt = merge_options(opt, varargin)
    for i = 1:2:numel(varargin)
        if isfield(opt, varargin{i})
            opt.(varargin{i}) = varargin{i+1};
        end
    end
end