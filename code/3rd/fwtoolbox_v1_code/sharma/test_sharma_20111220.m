% Author:       Samir Sharma
% Created:      October 2011
% Last Update:  December 20th, 2011 

clear all; close all; clc;

%* BEGIN: Load data, error check *%
foldername = '../../fwtoolbox_v1_data/kellman_data/';
fn = dir([foldername '*.mat']);
file_index = ceil(rand*length(fn));
disp([foldername fn(file_index).name]);
load([foldername fn(file_index).name]);
try
    data = imDataParams;
end
%* END: Load data, error check *%


%* BEGIN: Set algorithm parameters *
algoParams.species(1).name = 'water';
algoParams.species(1).ppm = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).ppm = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

% algorithm specific parameters
algoParams.stepsize = 0.75;     % support size scaling 
algoParams.min_win_size = 16;   % minimum 1D support size for B-spline
algoParams.MaxIter = 10;        % maximum iterations per scale
%* END: Set algorithm parameters *


%* BEGIN: Iterative decomposition *
[outParams] = fw_i2cm0i_3pluspoint_sharma(data,algoParams);
%* END: Iterative decomposition *
