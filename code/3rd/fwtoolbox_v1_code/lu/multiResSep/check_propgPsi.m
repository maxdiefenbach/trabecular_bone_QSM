function psiHat = check_propgPsi(psiHat, pp, contourInd, mask, swWidth)
[ht, wd] = size(psiHat);

% for each field map value, check its neighbors
if (exist('swWidth')==0)
swWidth = ceil(min(ht, wd)/10);
end
swWidth = min(10, swWidth);

for ctr=1:size(contourInd, 1),
    k = contourInd(ctr, 1); l = contourInd(ctr, 2);
    
    if (mask(k, l))
        
        swVrange = max(1,k-swWidth):min(ht, k+swWidth);
        swHrange = max(1,l-swWidth):min(wd, l+swWidth);
        
        localPsiBlk = psiHat(swVrange, swHrange);
        localMaskBlk= mask(swVrange, swHrange);
        localPsiBlk = localPsiBlk(localMaskBlk);
        localPsiBlk = localPsiBlk(:);
        localPsiMedian = median(localPsiBlk);
        
        currPsi = psiHat(k, l);
        localPsiDiff = max(abs(localPsiBlk-currPsi));
        if localPsiDiff>50,
            pset = pp(k, l, :); pset = pset(:);
            [diffPsi, minInd] = min(abs(pset - localPsiMedian));
            if diffPsi<50,
                psiHat(k, l) = pset(minInd);
            else
                psiHat(k, l) = localPsiMedian;
            end
        end
        
    end
end

return;