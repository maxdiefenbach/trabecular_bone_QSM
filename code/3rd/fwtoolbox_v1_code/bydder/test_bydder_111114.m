%% test_bydder_111110
%%
%% Test Mark Bydder's fat-water algorithms (water/fat modeling and phase-constrained water/fat separation) on Peter Kellman's data
%%
%% Author: Diego Hernando
%% Date created: August 18, 2011
%% Date last modified: November 14, 2011

% Add to matlab path
BASEPATH = '../hernando/';
addpath([BASEPATH 'common/']);
addpath([BASEPATH 'graphcut/']);
addpath([BASEPATH 'descent/']);
addpath([BASEPATH 'mixed_fitting/']);
addpath([BASEPATH 'create_synthetic/']);
mymachine = computer;
% LINUX
if strcmp(mymachine(1:3),'GLN')
    if strcmp(mymachine,'GLNX86') 
      addpath([BASEPATH 'matlab_bgl-4.0.1/']);
      disp ('adding bgl-4.0.1');
    else
      addpath([BASEPATH 'matlab_bgl-3.0-beta/']);
      disp ('adding bgl-3.0-beta');
    end
end

% PCWIN
if strcmp(mymachine(1:5),'PCWIN')
  addpath([BASEPATH 'matlab_bgl-wintemp/']);
  disp ('adding bgl wintemp');
end

% MAC
if strcmp(mymachine(1:3),'MAC')
  addpath([BASEPATH 'matlab_bgl_macosx64/']);
  disp ('adding bgl_mac');
end


%% Load some data
foldername = [BASEPATH '../../fwtoolbox_v1_data/kellman_data/'];
fn = dir([foldername '*.mat']);
file_index = ceil(rand*length(fn));
load([foldername fn(file_index).name]);
imDataParams = data;
imDataParams.images = double(data.images);
imDataParams.FieldStrength = 1.5;
%imDataParams.PrecessionIsClockwise = -1;

%% Set recon parameters
% General parameters
numberFatDoubleBonds = 3;
waterFrequencyPPM = 4.7;
algoParams.species = getFatWaterSpectrum_massScale_fatByDoubleBonds(numberFatDoubleBonds,waterFrequencyPPM);

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 100]; % Range of R2* values
algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD = 0.01;

%% Recon -- graph cut 
%% (Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)
tic
  outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );
toc

[outParamsCommonPhase] = decomposeGivenFieldMapAndDampings_commonPhase( imDataParams,algoParams,outParams.fieldmap,outParams.r2starmap,outParams.r2starmap );