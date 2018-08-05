%% test_MultiSeedRegionGrowing
%%
%% Author: Johan Berglund
%% Date created: November 16, 2011
%% Date last modified: November 18, 2011

% Add to matlab path
BASEPATH = './';
addpath([BASEPATH 'Common/']);
addpath([BASEPATH 'MultiSeedRegionGrowing/']);

%% Load some data
%foldername = [BASEPATH '../../New_Kellman_data/'];
foldername = [BASEPATH '../../fwtoolbox_v1_data/kellman_data/'];
fn = dir([foldername '*.mat']);
file_index = ceil(rand*length(fn));
load([foldername fn(file_index).name]);
imDataParams = data;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.70;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
algoParams.species(2).relAmps = [88 642 58 62 58 6 39 10 37];

% Algorithm-specific parameters
algoParams.c1 = 0.75; % Magnitude weight threshold for seed points
algoParams.c2 = 0.25; % Threshold on |log(W/F)| for seed points

%% Recon -- multi-seed region growing scheme for three-point data
%% (Berglund J, Johansson L, Ahlström H, Kullberg J. Three-point Dixon method enables whole-body 
%%  water and fat imaging of obese subjects. Magn Reson Med. 2010 Jun;63(6):1659-1668.)
tic
  outParams = fw_i3cm0i_3point_berglund( imDataParams, algoParams );
toc