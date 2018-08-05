%%
% Author: Harry Hu
% Date of last update: December 14, 2011
% SEE ALSO: README.PPT file for flowchart

%% 
% This folder contains a main program called "fw_ic2m0c_3pluspoint.m" and 
% is the implementation of region growing IDEAL for fat-water separation. 
% Program is based on the papers:
%   Yu, et al. MRM 2005 54:1032-1039 
%   Reeder, et al. MRM 2004 51:35-45

% Code does not perform R2* estimation
% Code can handle any # of echoes, preferably >= 3 for fw separation
% stability.

%%
% It's highly recommended that one reads the above papers to understand
% the syntax of the program.  I've tried to duplicate the algorithm as much
% as possible from Huanzhou Yu's paper.

%%
% test_fw_i2cm0c_3pluspoint_RGsinglecoil.m
%
% The code currently is written to reconstruct one coil per implementation.
% Thus if input contains more than one coil, program will ask which coil to
% reconstruct. 

% test_fw_i2cm0c_3pluspoint_RGmulticoil.m
%
% Alternatively one can implement coil combination with Diego Hernando's
% coilCombine function.

% Code can handle multiple slices.
%%
% Code written for teaching / demonstration purposes.
% Because of this, there are lots of interactive features.
% I use it as a teaching tool, and often ask students as a project to make 
% it run faster and more efficient.
% Not meant for large 3D data sets due to slow iterations.

%%
% The code illustrates some of the ad hoc nature of RG algorithm, 
% which requires a moderate amount of user adjustments to parameters.  
% These include: 
%     # of iterations
%     size of region growing kernel
%     how much downsample the data for initial RG loop

%%
% Code also implements voxel independent (VI) IDEAL (no RG), but note that 
% VI IDEAL can lead to severe water-fat swaps in complicated data sets.

%%
% The algorithm performs fairly well in data sets where the anatomy is 
% continuous.  However, region growing field map estimation fails in 
% anatomies such as the heart (data provided by Peter Kellman), due to 
% large air (no signal) areas (e.g. lungs), and lead to signal swaps.

%%
% Directory tree:
% ../ fw_i2cm0c_3plupoint.m
    % Should modify algoParams structure as needed based on region growing
    % Once implemented, code will add path for ../Support Functions and ../Data
% ../Support Functions
    % Ancillary support functions are all kept here
    % NO NEED TO MODIFY
% ../Data
    % Input .mat files should be kept here
    
%%
% The code is relatively simple Matlab syntax.  No Mex files.  It runs
% successfully on

% -- Linux, 64 and 32 bits (R2010b) 
% -- Mac OS X, 64 bit (2010b)
    %(v7.11.0.584,	64bit “maci64” R2010b)
% -- Windows XP and 7, 64 and 32 bits (R2010b 
    %(v7.10.0.499, 64bit R2010a)  
    %(v7.11.0.584, 32bit R2010b)
    
%%
% INPUT: matfile of data structure
% data.images [Nx Ny Nz Ncoils=1 NTE]
    % data.TE = TE values in (seconds)
    % data.FieldStrength in (Tesla)
    % data.PrecessionisClockwise (0 or 1)
% OUTPUT: matfile of output structure
    % outParams.water,  also outParams.species(1).amps
    % outParams.fat,    also outParams.species(2).amps
    % outParams.fieldmap (Hz)
    
%%
    
    