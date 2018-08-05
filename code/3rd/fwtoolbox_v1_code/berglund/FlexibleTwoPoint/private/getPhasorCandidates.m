%% Function name: getPhasorCandidates
%%
%% Description: Returns two candidates for the field map phasor, 
%%              following Eqs. 4-6 in:
%%
%% Berglund J, Ahlström H, Johansson L, Kullberg J. Two-point Dixon method
%% with flexible echo times. Magn Reson Med. 2011 Apr;65(4):994-1004.
%%
%% Input:
%%   - S: complex dual-echo image, size[nx,ny,nz,2]
%%   - a: complex linear model parameters, size[2]
%%
%% Output:
%%   - bA: complex field map phasor image 1st candidate, size[nx,ny,nz]
%%   - bB: complex field map phasor image 2nd candidate, size[nx,ny,nz]
%%
%% Author: Johan Berglund
%% Date created: November 17, 2011
%% Date last modified: November 18, 2011

function [bA bB] = getPhasorCandidates(S,a)
    %% Eq. 5:
    c1 = abs(S(:,:,:,1)).^2*(1-real(a(2)))-abs(S(:,:,:,2)).^2*(1-real(a(1)));
    c2 = abs(S(:,:,:,1)).^2*(abs(a(2))^2-real(a(2)))-abs(S(:,:,:,2)).^2*(abs(a(1))^2-real(a(1)));
    c3 = abs(S(:,:,:,1)).^2.*abs(S(:,:,:,2)).^2*abs(a(1)-a(2))^2-(imag(a(1))*abs(S(:,:,:,2)).^2-imag(a(2))*abs(S(:,:,:,1)).^2).^2;

    %% Eq. 4:
    QA = (c1+sqrt(c3))./(c1+c2);
    QB = (c1-sqrt(c3))./(c1+c2);

    %% Eq. 6:
    bA = S(:,:,:,2).*(1+QA*(a(1)-1))./S(:,:,:,1).*(1+QA*(a(2)-1));
    bB = S(:,:,:,2).*(1+QB*(a(1)-1))./S(:,:,:,1).*(1+QB*(a(2)-1));

    %% normalize results:
    bA = bA./abs(bA);
    bB = bB./abs(bB);

    %% make dummy assignment to any undefined voxels
    bA(isnan(bA))=1;
    bB(isnan(bB))=-1;
end