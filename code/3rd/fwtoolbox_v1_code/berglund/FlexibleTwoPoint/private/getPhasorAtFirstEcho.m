%% Function name: getPhasorAtFirstEcho
%%
%% Description: Returns the water phasor at the first echo, given the field
%%              map phasor, following Eq. 7 in:
%%
%% Berglund J, Ahlström H, Johansson L, Kullberg J. Two-point Dixon method
%% with flexible echo times. Magn Reson Med. 2011 Apr;65(4):994-1004.
%%
%% Input:
%%   - S: complex dual-echo image, size[nx,ny,nz,2]
%%   - a: complex linear model parameters, size[2]
%%   - b: complex field map phasor image, size[nx,ny,nz]
%%
%% Output:
%%   - b0: complex phasor image, size[nx,ny,nz]
%%
%% Author: Johan Berglund
%% Date created: November 17, 2011
%% Date last modified: November 18, 2011

function b0 = getPhasorAtFirstEcho(S,a,b)
    %% Eq. 7:
    b0 = (S(:,:,:,1)*(1-a(2))-S(:,:,:,2)*(1-a(1))./b)/(a(1)-a(2));
    
    %% Normalize result:
    b0 = b0./abs(b0);
end