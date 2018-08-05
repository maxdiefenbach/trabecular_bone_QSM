function [wat, fat, psiHat, resd] = multiResSep(TE, sigSet, df, dfAmp, vlevel)
% Inputs:
%        TE: msec (three echo times)
%        sigSet = cat(3, im1, im2, im3);
%        fieldStrength: in Tesla
%        df: fat frequencies (Hz)
%        dfAmp: fat relative amplitudes
%        vlevel: verbose level
%
%
% Outputs:
%        wat: water image
%        fat: fat image
%        psiHat: field map
%        resd: estimation map
%
% Wenmiao Lu 
% Stanford, May 2006
%

if (exist('vlevel')==0) vlevel = 1; end % default verbose level
if vlevel>0
disp('-------------------------------------------------');
disp('  Multi-resolution water and fat separation.     ');
disp('-------------------------------------------------');
end
%TE = TE(:)./1000;  % convert msec to sec.
[ht, wd, nechoes] = size(sigSet);
if (length(TE)~=nechoes),
    error('Echo times mismatch with data size.'); 
end

%% prepare the pyramids
% the coarsest resolution
minRes = 32; 
% number of levels
D = floor(log2(max(ht, wd)/minRes))+1;  
% initialize the pyramids
[imPyrm, watPyrm, fatPyrm, psiPyrm, resdPyrm] = ...
    init_pyramid(sigSet, D, ht, wd);

%% loop through pyramid levels
d = D;
while (D>0)
    if (vlevel>0)
       disp(['Level ' num2str(d)]);
    end
    
    if d==D, 
        % the coarsest level
        psiPyrm{d} = est_fieldMap(TE, imPyrm{d}, df, dfAmp , vlevel);          
    else
        % propagate field map to next fine level
        if d>0,
            psiPyrm{d} = propg_psi(psiPyrm{d+1}, TE, imPyrm{d}, df, dfAmp , vlevel);                
        else
            psiPyrm{d} = imresize(psiPyrm{d+1}, 2, 'bil');
        end
    end
    [watPyrm{d}, fatPyrm{d}, resdPyrm{d}] = ...
        lsSep_wf(TE, imPyrm{d}, psiPyrm{d}, df, dfAmp );
    
    if d==1,
        break; % the finest level done.
    else
        d = d - 1;
    end
end

%% final separation
% psiHat = psiPyrm{1}.*imMagMask(sigSet);
psiHat = psiPyrm{1};
ks = 5; psiHat = conv2(psiHat, ones(ks)/ks^2, 'same');
[wat, fat, resd] = lsSep_wf(TE, imPyrm{1}, psiHat, df, dfAmp );

if vlevel>0
figure; imshow(abs([wat fat]), []); title('Est. Water and Fat images');
figure; imshow(psiHat, []); title('Est. Field Map');
figure; imshow(resd, []); title('Estimation Error Map');
end

if nargout==1,
    wat = sum(sum(resd)); % only output the residue to reflect the performance.
end

return;
