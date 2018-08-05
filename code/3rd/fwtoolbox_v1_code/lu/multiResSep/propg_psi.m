function pp1 = propg_psi(pp2, TE, S, df, dfAmp , vlevel)
% 
% propagate the field from coarse to fine
%
% one coarse pixel corresponds to 4 fine pixels
%
% check the neighboring field map values, if there is a jump between
% neighboring values, use the nearest neighbor as initial starting value
%

mask = imMagMask(S(:,:,1), S(:,:,2), S(:,:,3));
[ht, wd] = size(mask);
% initial interpolation guess
[ht_pp2, wd_pp2] = size(pp2);
if (ht~=ht_pp2 | wd~=wd_pp2)
    pp2 = imresize(pp2, [ht wd], 'bil');
end
    
% water, fat, field map corrected images, residue
[w, f, resd2] = lsSep_wf(TE, S, pp2, df, dfAmp );
% compensate S with field map 
Sc = [];
for k=1:length(TE),
    Sc = cat(3, Sc, reshape(S(:,:,k).*exp(-i*2*pi*pp2.*TE(k)), [ht wd]));
end

A = calc_mixMat(TE, df, dfAmp );
deltaTE = abs(TE(2)-TE(1));
psiRange = 1/(4*deltaTE); % reduced search interval

pp1 = pp2; 
% refine pp2 to pp1 based on resd
resdThredHigh = .5;
resdThredLow = .1;
magMask = mean(abs(S), 3);

startInd = get_startPixInd(pp1, mask);
startIndy = mod(startInd, ht); if startIndy==0, startIndy=ht; end
startIndx = floor(startInd/ht); if startIndx==0, startIndx=1; end

load CONTOUR_IND; % load spiral contour indices
contourInd(:, 1) = contourInd(:, 1) + startIndy;
contourInd(:, 2) = contourInd(:, 2) + startIndx;
% get rid of outbounded regions
outBoundInd = find(contourInd(:, 1)<1 | contourInd(:, 1)>ht);
contourInd(outBoundInd, :) = [];
outBoundInd = find(contourInd(:, 2)<1 | contourInd(:, 2)>wd);
contourInd(outBoundInd, :) = [];

% neighboring field map values check and smoothing
pp1Init = pp1; pp1 = pp1.*0;
while (max(max(pp1-pp1Init)))>1
    pp1 = pp1Init;
    pp1Init = check_propgPsi(pp1, repmat(pp1, [1 1 3]), contourInd, mask);
end

for ctr = 1:size(contourInd, 1), 
          
    k = contourInd(ctr, 1); l = contourInd(ctr, 2);
    
    if (resd2(k, l) > resdThredLow*magMask(k,l) && mask(k, l)),
    
        swWidth = 5;
        swVrange = max(1,k-swWidth):min(ht, k+swWidth);
        swHrange = max(1,l-swWidth):min(wd, l+swWidth);
        localPsiBlk = pp1(swVrange, swHrange);
        localMaskBlk = mask(swVrange, swHrange);
        localPsiBlk = localPsiBlk(find(localMaskBlk~=0));
        localPsiBlk = localPsiBlk(:);
        localPsiMedian = mean(localPsiBlk);                                 
            
        tol = 1e-3; % threshold for golden section search
        if ( resd2(k, l)<resdThredHigh*magMask(k,l) & ...
                max(abs(pp1(k, l)-localPsiBlk))<15 ),
            
                % Sc contains field map compensated data 
                s = Sc(k, l, :); s = s(:);                             
                
                dpsi = golden('rfunc', ...
                    -1/2*psiRange, 0, 1/2*psiRange, ...
                    tol,s,A,TE, vlevel);
                
                % refine previous estimate
                pp1(k, l) = pp1(k, l) + dpsi; 
                
            else

                % S contains uncompensated data
                s = S(k, l, :); s = s(:);                               
        
                dpsi = golden('rfunc', ...
                    localPsiMedian-1/2*psiRange, ...
                    localPsiMedian, ...
                    localPsiMedian+1/2*psiRange, ...
                    tol,s,A,TE, vlevel);

                pp1(k, l) = dpsi; 
               
        end

    end
end

pp1Init = pp1; pp1 = pp1.*0;
while (max(max(pp1-pp1Init)))>1
    pp1 = pp1Init;
    pp1Init = check_propgPsi(pp1, repmat(pp1, [1 1 3]), contourInd, mask);
end
% pp1 = pp1.*mask;

return;