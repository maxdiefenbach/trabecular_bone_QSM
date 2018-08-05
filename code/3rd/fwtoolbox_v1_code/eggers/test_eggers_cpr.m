%% test_eggers_cpr
%%
%% Test conjugate phase reconstruction on (Holger Eggers') data
%%
%% Author: Holger Eggers
%% Date created: January 11, 2012
%% Date last modified: January 11, 2012

% Set path
BASEPATH = './';

addpath([BASEPATH '../hernando/common/']);
addpath([BASEPATH '../hernando/graphcut/']);
addpath([BASEPATH '../hernando/descent/']);
addpath([BASEPATH '../hernando/mixed_fitting/']);
addpath([BASEPATH '../hernando/create_synthetic/']);
mymachine = computer;
% LINUX
if strcmp(mymachine(1:3),'GLN')
    if strcmp(mymachine,'GLNX86') 
      addpath([BASEPATH '../hernando/matlab_bgl-4.0.1/']);
      disp ('adding bgl-4.0.1');
    else
      addpath([BASEPATH '../hernando/matlab_bgl-3.0-beta/']);
      disp ('adding bgl-3.0-beta');
    end
end

% PCWIN
if strcmp(mymachine(1:5),'PCWIN')
  addpath([BASEPATH '../hernando/matlab_bgl/']);
  disp ('adding bgl wintemp');
end

% MAC
if strcmp(mymachine(1:3),'MAC')
  addpath([BASEPATH '../hernando/matlab_bgl_macosx64/']);
  disp ('adding bgl_mac');
end


% Load input image data
foldername = [BASEPATH '../../fwtoolbox_v1_data/eggers_data_spiral/'];
load([foldername 'HEdata_spiral.mat']);

imDataParams        = data;
imDataParams.images = double(data.images);
%imDataParams.PrecessionIsClockwise = 1;

% Set parameters for water-fat separation (copied from Diego Hernando's test program)
algoParams.species(1).name       = 'water';
algoParams.species(1).frequency  = 0;
algoParams.species(1).relAmps    = 1;
algoParams.species(2).name       = 'fat';
algoParams.species(2).frequency  = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps    = [0.087 0.693 0.128 0.004 0.039 0.048];

algoParams.size_clique           = 1;          % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star          = [0 0];      % Range of R2* values
algoParams.NUM_R2STARS           = 1;          % Number of R2* values for quantization
algoParams.range_fm              = [-400 400]; % Range of field map values
algoParams.NUM_FMS               = 301;        % Number of field map values to discretize
algoParams.NUM_ITERS             = 40;         % Number of graph cut iterations
algoParams.SUBSAMPLE             = 2;          % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT                 = 0;          % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER            = 2;          % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda                = 0.05;       % Regularization parameter
algoParams.LMAP_EXTRA            = 0.05;       % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD                        = 0.01;

% Run water-fat separation (copied from Diego Hernando's test program)
tic
  outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );
toc

% Display intermediate results
figure(1);
colormap('gray');
imagesc(abs(outParams.species(1).amps));
figure(2);
colormap('gray');
imagesc(abs(outParams.species(2).amps));
figure(3);
imagesc(outParams.fieldmap);
colorbar;

% Set parameters for field inhomogeneity correction
algoParamsCpr.nv      =  15;         % Number of spiral interleaves
algoParamsCpr.ns      =  6434;       % Number of spiral samples per interleaf

algoParamsCpr.alpha   =  1.25;       % Gridding oversampling factor

algoParamsCpr.ks      =  4;          % Gridding kernel size
algoParamsCpr.gti     =  1000;       % Gridding kernel table increment

algoParamsCpr.tau     =  2.325e-6;   % Sampling interval [s] 

algoParamsCpr.shutter =  1;          % Circular image space shutter

% Compute trajectory
is                    =  size( outParams.species(1).amps, 1 );

lambda                =  3.0;        % Shape parameter
delay                 =  0.0;

trajectory =  compute_constant_density_spiral_trajectory( algoParamsCpr.nv, algoParamsCpr.ns, is, lambda, delay );

% Compute sampling density compensation weights
density_compensation =  compute_constant_density_spiral_density( algoParamsCpr.nv, algoParamsCpr.ns, 2 * ceil( algoParamsCpr.alpha * is / 2 ), is, lambda, delay, trajectory );

trajectory =  is * trajectory + is / 2;

algoParamsCpr.trajectory           =  trajectory; 
algoParamsCpr.density_compensation =  density_compensation;

gamma            = 42.6e6;
chemical_shift   = -3.4e-6; % Fat modeled as single peak here to keep it simple
frequency_offset = gamma * chemical_shift * imDataParams.FieldStrength;

outParams.species(1).fieldmap = - outParams.fieldmap;
outParams.species(2).fieldmap = - outParams.fieldmap - frequency_offset;

% Run field inhomogeneity correction
% (Eggers H, Knopp T, Potts D. Field inhomogeneity correction based on gridding reconstruction for 
%  magnetic resonance imaging. IEEE Trans Med Imaging 2007; 26:374-384.
tic
  outParamsCpr = cpr( outParams, algoParamsCpr );
toc

% Display results
figure(4);
colormap('gray');
imagesc(abs(outParamsCpr.species(1).amps));
figure(5);
colormap('gray');
imagesc(abs(outParamsCpr.species(2).amps));
