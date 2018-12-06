function [xMin, xMax] = getMinMax(x, nn)

    % Prepare sparse matrix vectors
    num = numel(nn);
    ii  = [rldecode((1:num)', nn, 1);
           rldecode((1:num)', max(nn) - nn, 1)];
    jj  = [mcolon(1, nn)';
           mcolon(nn+1, max(nn))'];
    nanVec = nan(sum(max(nn) - nn), 1);
    
    % Compute minimum and maximum cell coordinates
    [xMin, xMax] = deal(zeros(num, size(x,2)));
    for dNo = 1:size(x,2)
        v = [x(:,dNo); nanVec];
        xMat = sparse(ii, jj, v, num, max(nn));
        xMin(:,dNo) = min(xMat, [], 2);
        xMax(:,dNo) = max(xMat, [], 2);
    end
    
end