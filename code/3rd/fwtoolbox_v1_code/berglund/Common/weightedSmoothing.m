%% Function name: weightedSmoothing
%%
%% Description: Weighted smoothing of phasor image by 3D average filter
%%
%% Input:
%%   - phasor: complex phasor image, size[nx,ny,nz]
%%   - weight: weight image, size[nx,ny,nz]
%%   - voxelSize: voxel dimension in mm, size[3]
%%   - kernel: filter kernel size in mm, size[1]
%%
%% Output:
%%   - res: smooth complex phasor image, size[nx,ny,nz]
%%
%% Author: Johan Berglund
%% Date created: November 17, 2011
%% Date last modified: November 18, 2011

function res = weightedSmoothing(phasor,weight,voxelSize,kernel)
    %% translate kernel(mm) to span(voxels) in all directions and assert odd spans
    span = round((kernel./voxelSize-1)/2)*2+1;
    
    %% weighting of phasor image
    res = phasor.*weight;
    
    %% 3D average filter
    res=filter(ones(1,span(1))/span(1),1,res,[],1); % x direction
    res=filter(ones(1,span(2))/span(2),1,res,[],2); % y direction
    res=filter(ones(1,span(3))/span(3),1,res,[],3); % z direction
    
    %% normalize result
    res = res./abs(res);   
end