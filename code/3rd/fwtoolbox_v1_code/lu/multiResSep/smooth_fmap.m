function psiHat = smooth_fmap(pInitMap, pp, psiRange, mask, vlevel)

[ht, wd] = size(mask);

%% use the global median value to ensure the smoothness of the field map
% psiHatSet = repmat(zeros(size(pInitMap)), [1 1 3]);
% for ii = 1:3,    
%     
%     deltaPsiMedian = inf;
%     psiMedian = inf;
%     currPsiHat = pInitMap;
%     while deltaPsiMedian>.1
%         
%         oldPsiMedian = psiMedian;
%         
%         pm = currPsiHat(mask); pm = sort(pm(:));
%         if ii==1, 
%             psiMedian = median(pm(1:round(length(pm)/2)));
%         elseif ii==2,
%             psiMedian = median(pm);
%         else
%             psiMedian = median(pm(1+round(length(pm)/2):end));
%         end
%         
%         for k=1:ht,
%             for l=1:wd,        
%                 pset = pp(k, l, :); pset = pset(:);
%                 [diffPsi, minInd] = min(abs(pset - psiMedian));
%                 currPsiHat(k, l) = pset(minInd);
%             end
%         end        
%         
%         deltaPsiMedian = abs(psiMedian-oldPsiMedian);
%     end
%     psiHatSet(:,:,ii) = currPsiHat;
% end
psiHatSet = pp;

%% local psi fitting
% psiHat in the previous step gives the indication of reliable psi
% estimation, starting from that, we perform local median shift
% for each map, measure its quality, and pick the best
psiHat = bestInitMap(psiHatSet, mask);
% identify the reliable psi estimation from psiHat
% [startIndx, startIndy] = get_startCentroidXY(mask);
[startInd, startPsi] = get_startPixInd(psiHat, mask);
startIndy = mod(startInd, ht); if startIndy==0, startIndy=ht; end
startIndx = ceil(startInd/ht); if startIndx==0, startIndx=1; end
psiHat(startIndy, startIndx) = startPsi;
if vlevel>2,
r = mask; g = mask; b = mask;
g(startIndy, startIndx) = 0; 
b(startIndy, startIndx) = 0;
dum = uint8(cat(3, r, g, b));
figure; imshow(dum.*256, []); title('starting location');
end

%% progressive field map growth
diffPsiThresh = 100;
for k=1:max(ht, wd)
    
    % contour 1 is the starting pixel
    % first, get 4 corners for next contour
    topLeftCorner = [max(1, startIndy-k) max(1, startIndx-k)];
    topRightCorner= [max(1, startIndy-k) min(wd, startIndx+k)];
    botLeftCorner = [min(ht, startIndy+k) max(1, startIndx-k)];
    botRightCorner= [min(ht, startIndy+k) min(wd, startIndx+k)];
    
    % second, set the four corners
    % 1) topLeftCorner
    currIndy = topLeftCorner(1); currIndx = topLeftCorner(2);
    pset = pp(currIndy, currIndx, :); pset = pset(:);
    [diffPsi, minInd] = min(abs(pset - psiHat(currIndy+1, currIndx+1)));
    if diffPsi<diffPsiThresh,    
        psiHat(currIndy, currIndx) = pset(minInd);  
    else
        psiHat(currIndy, currIndx) = psiHat(currIndy+1, currIndx+1);
    end
    % 2) topRightCorner
    currIndy = topRightCorner(1); currIndx = topRightCorner(2);
    pset = pp(currIndy, currIndx, :); pset = pset(:);
    [diffPsi, minInd] = min(abs(pset - psiHat(currIndy+1, currIndx-1)));
    if diffPsi<diffPsiThresh,    
        psiHat(currIndy, currIndx) = pset(minInd);  
    else
        psiHat(currIndy, currIndx) = psiHat(currIndy+1, currIndx-1);
    end
    % 3) botLeftCorner
    currIndy = botLeftCorner(1); currIndx = botLeftCorner(2);
    pset = pp(currIndy, currIndx, :); pset = pset(:);
    [diffPsi, minInd] = min(abs(pset - psiHat(currIndy-1, currIndx+1)));
    if diffPsi<diffPsiThresh,    
        psiHat(currIndy, currIndx) = pset(minInd);  
    else
        psiHat(currIndy, currIndx) = psiHat(currIndy-1, currIndx+1);
    end
    % 4) botRightCorner
    currIndy = botRightCorner(1); currIndx = botRightCorner(2);
    pset = pp(currIndy, currIndx, :); pset = pset(:);
    [diffPsi, minInd] = min(abs(pset - psiHat(currIndy-1, currIndx-1)));
    if diffPsi<diffPsiThresh,    
        psiHat(currIndy, currIndx) = pset(minInd); 
    else
        psiHat(currIndy, currIndx) = psiHat(currIndy-1, currIndx-1);
    end
    
    % third, set the boundaries of each contour in turns
    % 1) top boundary
    currIndy = topLeftCorner(1);
    for l=topLeftCorner(2)+1:topRightCorner(2)-1,        
        currIndx = l;
        pset = pp(currIndy, currIndx, :); pset = pset(:);
%         [diffPsi, minInd] = min(abs(pset - psiHat(currIndy+1, currIndx)));        
        if k>1,
            dum = mean(psiHat(currIndy+1, l-1:l+1));
        else
            dum = psiHat(currIndy+1, currIndx);
        end
        [diffPsi, minInd] = min(abs(pset - dum));
        if diffPsi<diffPsiThresh,         
            psiHat(currIndy, currIndx) = pset(minInd);  
        else
            psiHat(currIndy, currIndx) = dum;
        end
    end
    % 2) left boundary
    currIndx = topLeftCorner(2);
    for l=topLeftCorner(1)+1:botLeftCorner(1)-1,
        currIndy = l;         
        pset = pp(currIndy, currIndx, :); pset = pset(:);
%         [diffPsi, minInd] = min(abs(pset - psiHat(currIndy, currIndx+1)));
        if k>1,        
            dum = mean(psiHat(l-1:l+1, currIndx+1));
        else
            dum = psiHat(currIndy, currIndx+1);
        end
        [diffPsi, minInd] = min(abs(pset - dum));
        if diffPsi<diffPsiThresh,         
            psiHat(currIndy, currIndx) = pset(minInd);  
        else
            psiHat(currIndy, currIndx) = dum;
        end
    end
    % 3) right boundary
    currIndx = topRightCorner(2);    
    for l=topRightCorner(1)+1:botRightCorner(1)-1,
        currIndy = l;         
        pset = pp(currIndy, currIndx, :); pset = pset(:);
%         [diffPsi, minInd] = min(abs(pset - psiHat(currIndy, currIndx-1)));
        if k>1,
            dum = mean(psiHat(l-1:l+1, currIndx-1));
        else
            dum = psiHat(currIndy, currIndx-1);
        end
        [diffPsi, minInd] = min(abs(pset - dum));
        if diffPsi<diffPsiThresh,
            psiHat(currIndy, currIndx) = pset(minInd);
        else
            psiHat(currIndy, currIndx) = dum;
        end
    end
    % 4) bottom boundary
    currIndy = botLeftCorner(1);
    for l=botLeftCorner(2)+1:botRightCorner(2)-1,
        currIndx = l;
        pset = pp(currIndy, currIndx, :); pset = pset(:);
%         [diffPsi, minInd] = min(abs(pset - psiHat(currIndy-1, currIndx)));
        if k>1,
            dum = mean(psiHat(currIndy-1, l-1:l+1));
        else
            dum = psiHat(currIndy-1, currIndx);
        end
        [diffPsi, minInd] = min(abs(pset - dum));
        if diffPsi<diffPsiThresh,         
            psiHat(currIndy, currIndx) = pset(minInd);  
        else
            psiHat(currIndy, currIndx) = dum;
        end
     end
    
%      imshow(psiHat, []); drawnow; 
end

%% check neighboring pixels to avoid abrupt changes
load CONTOUR_IND; % load spiral contour indices
contourInd(:, 1) = contourInd(:, 1) + startIndy;
contourInd(:, 2) = contourInd(:, 2) + startIndx;
% get rid of outbounded regions
outBoundInd = find(contourInd(:, 1)<1 | contourInd(:, 1)>ht);
contourInd(outBoundInd, :) = [];
outBoundInd = find(contourInd(:, 2)<1 | contourInd(:, 2)>wd);
contourInd(outBoundInd, :) = [];

pp = cat(3, psiHatSet-psiRange, psiHatSet, psiHatSet+psiRange);
psiHatInit = psiHat; psiHat = psiHat.*0;
while (max(abs(psiHatInit(:)-psiHat(:)))>1)
    psiHat = psiHatInit;
    psiHatInit = check_propgPsi(psiHat, pp, contourInd, mask, 1);    
end

return;