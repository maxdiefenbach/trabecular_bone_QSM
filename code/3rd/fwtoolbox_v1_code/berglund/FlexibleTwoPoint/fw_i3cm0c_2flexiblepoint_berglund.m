%% Function name: fw_i3cm0c_2flexiblepoint_berglund
%%
%% Description: Fat-water separation from two complex echoes with arbitrary
%%              choice of echo times, as described in:
%%
%% Berglund J, Ahlström H, Johansson L, Kullberg J. Two-point Dixon method
%% with flexible echo times. Magn Reson Med. 2011 Apr;65(4):994-1004.
%% 
%% Some properties:
%%   - Image-space
%%   - 2 species (water-fat)
%%   - Complex-fitting
%%   - Multi-peak fat (pre-calibrated)
%%   - No R2*
%%   - Common initial water/fat phase
%%   - Requires 2 flexible echoes
%%
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,nz,ncoils,2]
%%   - imDataParams.TE: echo times (in seconds)
%%   - imDataParams.FieldStrength: (in Tesla)
%%   - imDataParams.voxelSize: (mm x mm x mmm)
%%
%%   - algoParams.species(ii).name = name of species ii (string)
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude of each peak within species ii
%%
%%   Example
%%   - algoParams.species(1).name = 'water';
%%   - algoParams.species(1).frequency = 4.70;
%%   - algoParams.species(1).relAmps = 1;
%%   - algoParams.species(2).name = 'fat';
%%   - algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
%%   - algoParams.species(2).relAmps = [88 642 58 62 58 6 39 10 37];
%%
%% Output: structure outParams
%%   - outParams.water: estimated water image, size [nx,ny,nz]
%%   - outParams.fat: estimated fat image, size [nx,ny,nz]
%%   - outParams.fieldmap: field map (in Hz, size [nx,ny,nz])
%%   - outParams.phase: phase of water at first echi (in rad, size [nx,ny,nz])
%%
%%
%% Author: Johan Berglund
%% Date created: November 17, 2011
%% Date last modified: December 20, 2011

function outParams = fw_i3cm0c_2flexiblepoint_berglund( imDataParams, algoParams )

%% Check validity of params, and set default algorithm parameters if not provided
[validParams,imDataParams,algoParams] = checkParamsAndSetDefaults_FlexibleTwoPoint( imDataParams,algoParams );
if validParams==0
  disp(['Exiting -- data not processed']);
  outParams = [];
  return;
end

%% create linear model matrix
%gyromagnetic ratio for hydrogen [MHz/T]
gyro = 42.6;
%resonance frequency vector [in radians]:
if (imDataParams.PrecessionIsClockwise<=0) %if not clockwise, (most) fat frequencies will be negative
    omega_p = -2*pi*gyro*imDataParams.FieldStrength*(algoParams.species(1).frequency-algoParams.species(2).frequency);
else %if clockwise, (most) fat frequencies will be positive
    omega_p = +2*pi*gyro*imDataParams.FieldStrength*(algoParams.species(1).frequency-algoParams.species(2).frequency);
end
a=exp(complex(0,imDataParams.TE(1:2)'*omega_p))*algoParams.species(2).relAmps';
A=[1 real(a(1)); 0 imag(a(1)); 1 real(a(2)); 0 imag(a(2))];

%% DH*: Combine multiple coils (replaced dummyCoilCombine)
S = dummyCoilCombine3D( imDataParams.images(:,:,:,:,1:2) );
if size(imDataParams.images,4) > 1
  S = permute(coilCombine3D( imDataParams.images(:,:,:,:,1:2) ),[1 2 3 5 4]);
else
  S = permute( imDataParams.images(:,:,:,:,1:2),[1 2 3 5 4]);
end

%% get two field map phasor candidates in each voxel
[bA bB] = getPhasorCandidates(S,a);

%% choose one candidate in each voxel:
b = trivialChoiceBetweenPhasorCandidates(bA, bB);

%% weighted smoothing of field map phasor
% magnitude image for weighted smoothing
magn=abs(S(:,:,:,1))+abs(S(:,:,:,2));
% smoothing kernel size in mm
kernel = 10;
b = weightedSmoothing(b,magn,imDataParams.voxelSize,kernel);

%% get water phasor at first echo
b0 = getPhasorAtFirstEcho(S,a,b);

%% weighted smoothing of water phasor at first echo
b0 = weightedSmoothing(b0,magn,imDataParams.voxelSize,kernel);

%% remove "phase errors"
S(:,:,:,1)=S(:,:,:,1)./b0;
S(:,:,:,2)=S(:,:,:,2)./b0./b;

%% find least squares estimates for water and fat
[w,f]=SolveLS_2point(S,A);

%% Put results in outParams structure
outParams.species(1).name = 'water';
outParams.species(2).name = 'fat';
outParams.species(1).amps = abs(w);
outParams.species(2).amps = abs(f);
if (imDataParams.PrecessionIsClockwise<=0)
    outParams.fieldmap = +angle(b)/(2*pi*diff(imDataParams.TE(1:2)));
    outParams.phase = +angle(b0);
else
    outParams.fieldmap = -angle(b)/(2*pi*diff(imDataParams.TE(1:2)));
    outParams.phase = -angle(b0);
end