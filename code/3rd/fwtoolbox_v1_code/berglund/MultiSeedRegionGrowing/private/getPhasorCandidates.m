%% Function name: getPhasorCandidates
%%
%% Description: Returns two candidates for the field map phasor as in:
%%
%% Berglund J, Ahlström H, Johansson L, Kullberg J. Closed-form solution 
%% for the three-point Dixon method with advanced spectrum modeling. 
%% In: Proc. ISMRM 2011
%%
%% Input:
%%   - S: complex triple-echo image, size[nx,ny,nz,3]
%%   - A: complex linear model matrix, size[3,2]
%%
%% Output:
%%   - bA: complex field map phasor image 1st candidate, size[nx,ny,nz]
%%   - bB: complex field map phasor image 2nd candidate, size[nx,ny,nz]
%%
%% Author: Johan Berglund
%% Date created: November 16, 2011
%% Date last modified: November 18, 2011

function [bA bB] = getPhasorCandidates(S, A)
    c1 = A(1,1)*A(3,2)-A(3,1)*A(1,2);
    c2 = A(1,1)*A(2,2)-A(2,1)*A(1,2);
    c3 = A(2,1)*A(3,2)-A(3,1)*A(2,2);

    bA=(S(:,:,:,2)*c1+sqrt(S(:,:,:,2).*S(:,:,:,2)*c1*c1-4*S(:,:,:,1).*S(:,:,:,3)*c2*c3))./(2*S(:,:,:,1)*c3);
    bB=(S(:,:,:,2)*c1-sqrt(S(:,:,:,2).*S(:,:,:,2)*c1*c1-4*S(:,:,:,1).*S(:,:,:,3)*c2*c3))./(2*S(:,:,:,1)*c3);

    %% make dummy assignment to any undefined voxels
    bA(isnan(bA))=1;
    bB(isnan(bB))=-1;
end