%% Function name: SolveLS
%%
%% Description: Finds complex-valued least-squares estimates of water and fat,
%%              following Eq. 12 in:
%%
%% Berglund J, Johansson L, Ahlström H, Kullberg J. Three-point Dixon method enables whole-body 
%% water and fat imaging of obese subjects. Magn Reson Med. 2010, Jun;63(6):1659-1668.
%%
%% Input:
%%   - S: phase-corrected complex triple-echo image, size[nx,ny,nz,3]
%%   - A: complex-valued linear model matrix, size[3,2]
%%
%% Output:
%%   - W: complex-valued water image, size[nx,ny,nz]
%%   - F: complex-valued fat image, size[nx,ny,nz]
%%
%% Author: Johan Berglund
%% Date created: November 16, 2011
%% Date last modified: November 18, 2011

function [W,F] = SolveLS(S,A)
    %% Find pseudoinverse of A
    Ainv = pinv(A);

    [nx ny nz ~] = size(S);
    
    %% W and F according to Eq. 12
    W=reshape(reshape(S,[],3)*Ainv(1,:).',[nx ny nz]);
    F=reshape(reshape(S,[],3)*Ainv(2,:).',[nx ny nz]);
end