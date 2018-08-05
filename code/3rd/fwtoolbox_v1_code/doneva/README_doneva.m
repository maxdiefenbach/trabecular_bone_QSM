%% Author: Mariya Doneva
%% Date: December, 2011
%This folder contains functions for compressed sensing fat-water separation.
%
%% The main function is:
% 
% * fw_k2cs0i_3point_doneva_cswf.m
% Fat-water separation from undersampled data using compressed sensing
% 
% (Doneva M, Börnert P, Eggers H, Mertins A, Pauly J, Lustig M
%  Compressed sensing for chemical shift-based water-fat separation.
%  Magn Reson Med. 2010 Sep;64(6):1749-1759)
% 
%% Other important functions:
% Initialization:
% * methods/nonlin_cg.m
%   CS reconstruction for each TE. CS reconstruction method is based on
%   M.Lustig, D.Donoho and J.Pauly, "Sparse MRI: The application of compressed sensing for
%   rapid MR imaging" Magn Reson Med 2007 58: 1182-95 
% * methods/locatemins_outer1.m 
%   field map estimation computing possible field map candidates for each
%   pixel 
% * methods/rm2d2.m 
%   estimating a smooth field map by region merging
% 
% Integrated CS-WF separation:
% * methods/wfcs
%
%%
% * create_synthetic/createSynthetic_kSpace.m
% Creation of synthetic chemical shift-encoded datasets, based on a
% ground truth (fat, water, R2*, field map), model (fat peaks), and
% acquisition parameters (echo times, field strength)
% 
%%
% * Folder sampling_masks/ 
%   contains example sampling masks for three echo acquisition.
%  
% * Folder common/ 
%  contains Fourier and sparsity transforms
% 
% * Folder parameter_initialization/ 
% contains scripts for parameter initialization
%
% * Folder methods/
% contains computation routines used for field map estimation, compressed
% sensing reconstruction on individual images, and integrated compressed
% sensing fat water separation
%
% * Folder WaveletMEX/ contains mex functions for 2D wavelet transform
% using WaveLab 
% The full WaveLab package is available at:
% http://www-stat.stanford.edu/~wavelab/Wavelab_850
%   
% The top-level folder contains two test scripts: 
% 
% * testSynthetic_doneva.m (create synthetic dataset based on Shepp-Logan phantom and test the compressed sensing water fat separation)
% 
% * test_doneva.m (test compressed sensing water fat separation on datasets from Peter Kellman)
% 
% 
% 
% 
%  
% * The code has been tested on the following platform:
% 
%     -- Windows, 64bit (R2009a)

