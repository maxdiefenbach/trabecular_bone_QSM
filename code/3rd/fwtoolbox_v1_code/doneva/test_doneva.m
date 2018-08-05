%% test_doneva_28122011
%%
%% Test compressed sensing fat-water separation on (Peter Kellman's) data
%%
%% Author: Mariya Doneva
%% Date created: December 20, 2011
%% Date last modified: December 28, 2011

% Add to matlab path
BASEPATH = './';%Q:\Work\matlab\water_fat_toolbox1\';
addpath([BASEPATH 'create_synthetic/']);
addpath([BASEPATH 'parameter_initialization/']);
addpath([BASEPATH 'sampling_masks/']);
addpath([BASEPATH 'common/']);
addpath([BASEPATH 'methods/']);
addpath([BASEPATH 'waveletMEX/']);
addpath([BASEPATH]);



%% Load test data 
%% Datasets 3,4,7,12 and 13 have 3 echoes
foldername = ['../../fwtoolbox_v1_data/kellman_data/'];%Q:\Work\matlab\WFToolbox\kellman_data\'];
% Please select a dataset and a corresponding sampling mask by uncommenting
% the corresponding lines
%% Dataset 3
 filename = 'PKdata3.mat';  
 load mask_3p_256_192_2d

%% Dataset 4
% filename = 'PKdata4.mat';  
% load mask_3p_192_256_2d

%% Dataset 7
% filename = 'PKdata7.mat';  
% load mask_3p_256_256_2d
% load mask_3p_256_256_3d_f3_3

%% Dataset 12
% filename = 'PKdata12.mat';  
% load mask_3p_256_192_2d

%% Dataset 13
% filename = 'PKdata13.mat';  
% load mask_3p_192_256_2d

load([foldername filename]);



%% Simulate undersampled k-space data
TE = data.TE;



%% create sampling mask 
size_params = size(data.images);
sx = size_params(1);
sy = size_params(2);
N  = size_params(5);



%% set data structure

% use precomputed sampling mask for the simulation of undersampled data
mask = mask_3p;
kDataParams0.kSpaceLocations = [];
kDataParams0.kSpaceTimes = [];
for k = 1:N
[X,Y] = find(mask(:,:,k) >0);
kDataParams0.kSpaceLocations = [kDataParams0.kSpaceLocations; [X-sx/2-1 Y-sy/2-1]];
kDataParams0.kSpaceTimes = [kDataParams0.kSpaceTimes; TE(k)*ones(size(X))];
end

kDataParams = kDataParams0;
for kt = 1:N
    kSpaceData = fft2c(data.images(:,:,1,1,kt));
    init_images(:,:,kt) = ifft2c(mask(:,:,kt).*kSpaceData);
    locations_idx = find(kDataParams.kSpaceTimes==TE(kt));    
    idx = sub2ind([sx,sy],kDataParams.kSpaceLocations(locations_idx,1)+ sx/2 +1, kDataParams.kSpaceLocations(locations_idx,2)+ sy/2 +1 );
    temp = kSpaceData(idx);
    kDataParams.kSpaceValues(locations_idx) = temp;
end

kDataParams.FOV = [sx;sy;1];
kDataParams.FieldStrength = data.FieldStrength; % 120426: Field Strength specified in data


dTE = TE(2)-TE(1);
% Algorithm-specific parameters
algoParams.dTE = dTE;
algoParams.num_iter_rm = 8; % Maximal number of iterations for region merging
algoParams.num_outer_iter_wfcs = 5;  % Number of outer iterations for wf compressed sensing
algoParams.num_inner_iter_wfcs = 4;   % Number of inner iterations for wf compressed sensing
algoParams.Thresh = 0.08; 
algoParams.largeFM = 1;   % Consider few 1/dTE periods for the field map estimation

%% The default regularization parameters can be changed through the
% following parameters
% Regularization parameters for CS on indiviudual echo images
algoParams.lambda_tv_cs         = 0.04; % TV regularization parameter for
algoParams.lambda_wavelet_wfcs  = 0.02; % wavelets regularization parameter for wfcs
% Regularization parameters for CS-WF 
algoParams.lambda_tv_wfcs       = 0.04; % TV regularization parameter for wfcs
algoParams.lambda_wavelet_cs    = 0.02; % wavelets regularization parameter for cs




%% set reconstruction parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];




%% Recon -- CS water-fat
%% (Doneva M, Börnert P, Eggers H, Mertins A, Pauly J, Lustig M
%% Compressed sensing for chemical shift-based water-fat separation. Magn
%% Reson Med. 2010 Sep;64(6):1749-1759.)
timestamp = tic;
outParams = fw_k2cs0i_3point_doneva_cswf( kDataParams, algoParams );
total_time = toc(timestamp);
disp('Total time'); disp(total_time);


figure; imshow(abs(outParams.species(1).amps),[]); title('Water image');
figure; imshow(abs(outParams.species(2).amps),[]); title('Fat image');
figure; imshow((outParams.fieldmap),[]); title('Field map');


figure; imshow(abs([init_images(:,:,1) init_images(:,:,2)  init_images(:,:,3) ]),[]); title('Initial echo images');




