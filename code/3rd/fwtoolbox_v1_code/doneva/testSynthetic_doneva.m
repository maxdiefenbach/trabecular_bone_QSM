%% testSynthetic_doneva
%%
%% Test CS fat-water algorithms on synthetic data
%%
%% Author: Mariya Doneva
%% Date created: October 3, 2011
%% Date last modified: October 3, 2011


% Add to matlab path
BASEPATH = './';%Q:\Work\matlab\water_fat_toolbox1\';
addpath([BASEPATH 'create_synthetic/']);
addpath([BASEPATH 'parameter_initialization/']);
addpath([BASEPATH 'sampling_masks/']);
addpath([BASEPATH 'common/']);
addpath([BASEPATH 'methods/']);
addpath([BASEPATH 'waveletMEX/']);
addpath([BASEPATH]);


%% Create synthetic data and set acquisition parameters
sx = 256;sy=sx;
imtest = phantom(sx);
threshold = 0.8;
N = 3; TEinit = -0.4e-3; dTE = 1.6e-3;
TE = TEinit + [0:N-1]*dTE;

trueParams.species(1).amps = imtest.*(imtest<threshold); % Water
trueParams.species(2).amps = imtest.*(imtest>=threshold); % Fat

[X,Y] = meshgrid(linspace(-1,1,sx),linspace(-1,1,sy));
trueParams.fieldmap = 20*randn*ones(sx,sy) + 40*randn*X + 40*randn*Y + 100*randn*X.^2 + 100*randn*Y.^2  + 100*randn*X.*Y.^2 + 100*randn*X.^3  + 100*randn*Y.^3;
trueParams.r2starmap = 0*imtest;


% use precomputed sampling mask for the simulation mask_3p
load mask_3p_256_256_3d


kDataParams0.kSpaceLocations = [];
kDataParams0.kSpaceTimes = [];
for k = 1:N
[X,Y] = find(mask_3p(:,:,k) >0);
kDataParams0.kSpaceLocations = [kDataParams0.kSpaceLocations; [X-sx/2-1 Y-sy/2-1]];
kDataParams0.kSpaceTimes = [kDataParams0.kSpaceTimes; TE(k)*ones(size(X))];
end

kDataParams0.FieldStrength = 1.5;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];


% Simulate data
kDataParams = createSynthetic_kSpace( kDataParams0, algoParams, trueParams );

% Algorithm-specific parameters
algoParams.dTE = dTE;
algoParams.num_outer_iter_wfcs = 5;  % Number of outer iterations for wf compressed sensing
algoParams.num_inner_iter_wfcs = 4;   % Number of inner iterations for wf compressed sensing
algoParams.Thresh = 0.03;    
algoParams.largeFM = 1;   % Consieder few 1/dTE periods for the field map estimation


algoParams.lambda_tv_cs         = 0.1; % TV regularization parameter for cs
algoParams.lambda_wavelet_cs    = 0.0; % wavelets regularization parameter
algoParams.lambda_tv_wfcs       = 0.04; % TV regularization parameter for wfcs
algoParams.lambda_wavelet_wfcs  = 0.0; % wavelets regularization parameter



%% Recon -- CS water-fat
%% (Doneva M, Börnert P, Eggers H, Mertins A, Pauly J, Lustig M
%% Compressed sensing for chemical shift-based water-fat separation. Magn
%% Reson Med. 2010 Sep;64(6):1749-1759.)
tic
outParams = fw_k2cs0i_3point_doneva_cswf( kDataParams, algoParams );
toc


figure; imshow(abs(outParams.species(1).amps),[]); title('Water image');
figure; imshow(abs(outParams.species(2).amps),[]); title('Fat image');
figure; imshow((outParams.fieldmap),[]); title('Field map');




