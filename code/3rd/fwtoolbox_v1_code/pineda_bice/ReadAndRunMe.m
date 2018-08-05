% File Explaining the NSA Code for the Fat/Water Toolbox
% Author: Emily Bice 
% Email: emily.bice@gmail.com
% Advisor: Angel Pineda
% Email: pineda@fullerton.edu
% Updated 2012-01-19

% We assume that a set of input and output parameters are in the
% MATLAB workspace.  In our case, we used 'test_hernando_110818'.


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
fn = dir([foldername '*data4.mat']);
file_index = 1;%ceil(rand*length(fn));
disp([foldername fn(file_index).name]);
load([foldername fn(file_index).name]);
imDataParams = data;
imDataParams.images = double(data.images);
% $$$ imDataParams.FieldStrength = 1.5;
% $$$ imDataParams.PrecessionIsClockwise = -1;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 0]; % Range of R2* values
algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 0; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
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

% For the first part of the NSA analysis, we look at the NSA as a function
% of fat fraction.  We assume a two-peak model where the fat and water are
% aligned at the echo at pi/4 radians.

% First we generate the NSA normalized to be between 0 and N points for
% the field map and no R2* model:

FMK_NSA = fw_nsa_2pluspoint(imDataParams, algoParams);
% For testing, we have plots of the resulting NSA.
ha=axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
text(0.5,1,'Theoretical NSA for a range of fat fractions','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)

pause %pressing any key will make MATLAB continue
close

% By adding arguments we can compute the NSA for the model with unknown
% Field Map, given in the structure simParams:
simParams = [];
simParams.fieldmap = 110;

IDEAL_NSA = fw_nsa_2pluspoint(imDataParams, algoParams, simParams);
ha=axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
text(0.5,1,'Theoretical NSA for a range of fat fractions; \psi = 110','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)

pause
close 


% And for the model including a single R2*:

simParams.r2star = 50;

IDEAL_R2_NSA = fw_nsa_2pluspoint(imDataParams, algoParams, simParams);
ha=axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
text(0.5,1,'Theoretical NSA for a range of fat fractions; \psi = 110, R_2^* = 50','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)

pause
close

%We also have a sample code which carries out a Monte Carlo using least
%squares estimation to test the NSA.  Note that for testing purposes we
%used less fat/water ratios and a small number of realizations.  With these
%parameters, the Monte Carlo took 10 seconds in my laptop.

IDEAL_R2_NSA_MC = fw_mcsim_2pluspoint(imDataParams, algoParams, simParams, 11, 200, 100);
ha=axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
text(0.5,1,'Monte Carlo NSA for a range of fat fractions; \psi = 110, R_2^* = 50','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)

pause
close

%We also take the data from a reconstructed image to create a map of the
%theoretical NSA for the image.  This took about 25 seconds in my laptop.

NSA_MAP = fw_nsa_map(imDataParams, algoParams, outParams);
title('Reconstruction NSA map')

pause
close
