%% Function: checkParamsAndSetDefaults_cswf
%%
%% Description: Check validity of input parameters and set defaults for unspecified parameters
%%
%% Input:
%%   - kDataParams: TEs, k-space data and field strength
%%   - algoParams: algorithm parameters
%%
%% Output:
%%   - validParams: binary variable (0 if parameters are not valid for this algorithm)
%%   - algoParams2: "completed" algorithm parameter structure (after inserting defaults for unspecified parameters)
%%
%%
%% Author: Mariya Doneva
%% Date created:
%% Date last modified:
%%

function [validParams,algoParams2] = checkParamsAndSetDefaults_cswf( kDataParams,algoParams )

algoParams2 = algoParams;
validParams = 1;


%% Check validity of provided data and recon parameters
if kDataParams.FOV(3) > 1
    disp('ERROR: 2D recon -- expected FOV = [SX;SY;1]')
    validParams = 0;
end

if length(algoParams.species) > 2
    disp('ERROR: Too many species -- reconstruction supports 2 chemical species (water and fat)')
    validParams = 0;
end

TE = unique(kDataParams.kSpaceTimes);

if (length(TE)~=3)
    disp('ERROR: Reconstruction supports only 3 echo measurements')
    validParams = 0;
end

if((TE(2)-TE(1))-(TE(3)-TE(2))>eps)
    disp('ERROR: Equidistant echo timed are required')
    validParams = 0;
end

%%   algoParams.largeFM = 1;  % 0,1 flag to enable larger interval for the field map estimation
try
    algoParams2.largeFM = algoParams.largeFM;
catch
    algoParams2.largeFM = 0;
end

try
    algoParams2.Thresh = algoParams.Thresh;
catch
    algoParams2.Thresh = 0.01;
end




