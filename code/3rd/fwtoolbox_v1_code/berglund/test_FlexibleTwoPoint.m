%% test_FlexibleTwoPoint
%%
%% Author: Johan Berglund
%% Date created: November 17, 2011
%% Date last modified: December 20, 2011

% Add to matlab path
BASEPATH = './';
addpath([BASEPATH 'Common/']);
addpath([BASEPATH 'FlexibleTwoPoint/']);

%% Load some data
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

%% Recon -- two point Dixon method with flexible echo times
%% (Berglund J, Ahlström H, Johansson L, Kullberg J. Two-point Dixon method
%%  with flexible echo times. Magn Reson Med. 2011 Apr;65(4):994-1004.)
tic
  outParams = fw_i3cm0c_2flexiblepoint_berglund( imDataParams, algoParams );
toc