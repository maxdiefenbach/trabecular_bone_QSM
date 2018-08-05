function contourInd = get_contourInd(ht, wd, startIndy, startIndx)
%
% generate the indices of one contour, given
% [ht, wd] and 
% pp: field map 
% mask: foreground mask

% startInd = get_startPixInd(pp, mask);
% startIndy = mod(startInd, ht); if startIndy==0, startIndy=ht; end
% startIndx = floor(startInd/ht); if startIndx==0, startIndx=1; end

contourInd = [startIndy startIndx];

% four flags indicate if the contour hits respective boundaries
maxContourWid = max([ht, wd]);
k = 1;
dirFlag = 1;
while (k<maxContourWid)    
    dum = contourInd(end, :);
    
    if (dirFlag)
        
        for n=1:k,
        dum(2) = dum(2)+1; % right
        contourInd = cat(1, contourInd, dum);
        end
        
        for n=1:k,
        dum(1) = dum(1)-1; % up
        contourInd = cat(1, contourInd, dum);
        end
    
        dirFlag = 0;
               
    else

        for n=1:k,
        dum(2) = dum(2)-1; % left
        contourInd = cat(1, contourInd, dum);
        end
        
        for n=1:k,
        dum(1) = dum(1)+1; % down
        contourInd = cat(1, contourInd, dum);
        end
        
        dirFlag = 1;        
        
    end
    
    k = k+1;    
    
end

% outBoundInd = find(contourInd(:, 1)<1 | contourInd(:, 1)>ht);
% contourInd(outBoundInd, :) = [];
% outBoundInd = find(contourInd(:, 2)<1 | contourInd(:, 2)>wd);
% contourInd(outBoundInd, :) = [];
% 
%figure; hold on;
%for k=1:size(contourInd, 1)
%    plot(contourInd(k, 2), contourInd(k, 1), 'x'); drawnow;
%end
