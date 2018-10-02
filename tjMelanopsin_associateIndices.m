function theIdx = tjMelanopsin_associateIndices(dataTraceIdx, annotIdx)

% Iterate over the passed indices
for ii = 1:length(annotIdx)
    tmp = find(dataTraceIdx == annotIdx(ii));
    if ~isempty(tmp)
        theIdx(ii) = tmp(1);
    else
        % Deal with the case where we can't find the index in the data
        % trace
        theCandidates = [annotIdx(ii)-100:annotIdx(ii)+100];
        tmp1 = NaN*ones(size(theCandidates));
        for cc = 1:length(theCandidates)
            tmp0 = find(dataTraceIdx == theCandidates(cc));
            if ~isempty(tmp0)
                tmp1(cc) = tmp0(1);
            end
        end
        the0 = tmp1;
        the0(isnan(tmp1)) = round(interp1(theCandidates(~isnan(tmp1)), tmp1(~isnan(tmp1)), theCandidates(isnan(tmp1))));
        theIdx(ii) = the0(find(theCandidates == annotIdx(ii)));
    end
end
