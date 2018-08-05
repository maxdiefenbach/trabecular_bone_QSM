% Function: decomposeGivenFieldMapAndDampings_commonPhase
%
% Description: estimate water/fat images given the nonlinear parameters, imposing common phase for water and fat
% 
% Parameters:
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%%   - imDataParams.TEs: echo times (in seconds)
%%   - imDataParams.fieldStrength: (in Tesla)
%%
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0] 
%%      - algoParams.species(1).relAmps = [1]   
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%% 
%
%  - fieldmap: the estimated B0 field map
%  - r2starWater: the estimated water R2* map
%  - r2starFat: the estimated fat R2* map
%
% Returns: structure "outParams":
%%   - outParams.species(ii).name: name of the species (taken from algoParams)
%%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils] 
%%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny])
%%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%
% Author: Mark Bydder
% Date created: November 14, 2011
% Date last modified: November 14, 2011


function [outParams] = decomposeGivenFieldMapAndDampings_commonPhase( imDataParams,algoParams,fieldmap,r2starWater,r2starFat )

gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:))*(imDataParams.FieldStrength)];
ampWater = algoParams.species(1).relAmps;
relAmps = algoParams.species(2).relAmps;
images = imDataParams.images;
t = imDataParams.TE;

sx = size(images,1);
sy = size(images,2);
N = size(images,5);
C = size(images,4);

relAmps = reshape(relAmps,1,[]);


B1 = zeros(N,2);
B = zeros(N,2);
for n=1:N
  B1(n,:) = [ampWater*exp(j*2*pi*deltaF(1)*t(n)),sum(relAmps(:).*exp(j*2*pi*deltaF(2:end)*t(n)))];
end

amps = zeros(sx,sy,2,C);
for kx =1:sx
  for ky=1:sy
    s = reshape( squeeze(images(kx,ky,:,:,:)), [C N]).';

    B(:,1) = B1(:,1).*exp(j*2*pi*fieldmap(kx,ky)*t(:) - r2starWater(kx,ky)*t(:));
    B(:,2) = B1(:,2).*exp(j*2*pi*fieldmap(kx,ky)*t(:) - r2starFat(kx,ky)*t(:));

    amps(kx,ky,:,:) = B\s;
    amps(kx,ky,:,:) = cls4(B,s);
  end
end

% Put results in outParams structure
try
  outParams.species(1).name = algoParams.species(1).name;
  outParams.species(2).name = algoParams.species(2).name;
catch
  outParams.species(1).name = 'water';
  outParams.species(2).name = 'fat';
end  

outParams.species(1).amps = amps(:,:,1,:);
outParams.species(2).amps = amps(:,:,2,:);
outParams.r2starmap = r2starWater;
outParams.fieldmap = fieldmap;

