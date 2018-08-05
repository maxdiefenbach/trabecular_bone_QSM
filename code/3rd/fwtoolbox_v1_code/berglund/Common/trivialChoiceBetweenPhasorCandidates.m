%% Function name: trivialChoiceBetweenPhasorCandidates
%%
%% Description: Chooses the phasor with the smallest phase voxel-by-voxel
%%
%% Input:
%%   - phasor1, phasor2: complex phasor images, size[nx,ny,nz]
%%
%% Output:
%%   - phasor: complex phasor image, size[nx,ny,nz]
%%
%% Author: Johan Berglund
%% Date created: November 17, 2011
%% Date last modified: November 18, 2011

function phasor = trivialChoiceBetweenPhasorCandidates(phasor1, phasor2)
    phasor = phasor2;
    smallest=abs(angle(phasor1))<abs(angle(phasor2));
    phasor(smallest)=phasor1(smallest);
end
