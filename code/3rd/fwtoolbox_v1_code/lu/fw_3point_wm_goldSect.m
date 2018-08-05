%% Function name: fw_3point_wm_goldSect
%%
%% Description: Fat-water separation using golden section search
%%
%% 
%% Some properties:
%%   - Image-space
%%   - 2 species (water-fat) single-peak
%%   - Complex-fitting
%%   - Independent water/fat phase
%%   - Requires 3+ echoes at arbitrary echo times 
%%
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%%   - imDataParams.TE: echo times (in seconds)
%%   - imDataParams.FieldStrength: (in Tesla)
%%
%% Output: structure outParams
%%   - outParams.species(ii).name: name of the species (taken from algoParams)
%%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils] 
%%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%%
%%

function outParams = fw_3point_wm_goldSect( imDataParams, algoParams)


% DH* 111219: If precession is clockwise (positive fat frequency) simply conjugate data
if imDataParams.PrecessionIsClockwise <= 0 
  imDataParams.images = conj(imDataParams.images);
  imDataParams.PrecessionIsClockwise = 1;
end

% DH* Check multi-coil or multi-slice
% Get data dimensions
[sx,sy,sz,C,N] = size(imDataParams.images);
% DH* If more than one slice, pick central slice
if sz > 1
  disp('Multi-slice data: processing central slice');
  imDataParams.images = imDataParams.images(:,:,ceil(end/2),:,:);
end
% DH* If more than one channel, coil combine
if C > 1
  disp('Multi-coil data: coil-combining');
  imDataParams.images = coilCombine(imDataParams.images);
end
 
% DH* 111219: Get fat frequencies and relative amplitudes
if nargin < 2
  df = [-420 -318 94]*imDataParams.FieldStrength/3;
  dfAmp = [0.75 0.17 0.08]; % amplitudes of multiple fat peaks
else
  df = (algoParams.species(2).frequency - algoParams.species(1).frequency(1))*42.58*imDataParams.FieldStrength;
  dfAmp = algoParams.species(2).relAmps;
end

%addpath('./multiResSep');
[wat, fat, psiHat] = multiResSep(imDataParams.TE, imDataParams.images, df, dfAmp, 0);

outParams.species(1).name = 'water';
outParams.species(1).amps = (wat); %%DH* 111219: return complex images, in case phase is desired
outParams.species(2).name = 'fat';
outParams.species(2).amps = (fat); %%DH* 111219: return complex images, in case phase is desired
outParams.fieldmap = psiHat;
