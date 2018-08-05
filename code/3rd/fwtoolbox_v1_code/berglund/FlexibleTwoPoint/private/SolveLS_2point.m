%% Function name: SolveLS_2point
%%
%% Description: Finds real-valued least-squares estimates of water and fat,
%%              following Eq. 10 in:
%%
%% Berglund J, Ahlström H, Johansson L, Kullberg J. Two-point Dixon method
%% with flexible echo times. Magn Reson Med. 2011 Apr;65(4):994-1004.
%%
%% Input:
%%   - S: phase-corrected complex dual-echo image, size[nx,ny,nz,2]
%%   - A: real-valued linear model matrix, size[4,2]
%%
%% Output:
%%   - w: real-valued water image, size[nx,ny,nz]
%%   - f: real-valued fat image, size[nx,ny,nz]
%%
%% Author: Johan Berglund
%% Date created: November 17, 2011
%% Date last modified: November 18, 2011

function [w,f] = SolveLS_2point(S,A)
    %% Find pseudoinverse of A
    Ainv = pinv(A);                                                             

    %% Rearrange complex-valued S into real-valued S0
    [nx ny nz ~] = size(S);
    S0 = [reshape(real(S(:,:,:,1)),[],1) reshape(imag(S(:,:,:,1)),[],1) reshape(real(S(:,:,:,2)),[],1) reshape(imag(S(:,:,:,2)),[],1)];

    %% w and f according to Eq. 10
    w=reshape(S0*Ainv(1,:).',[nx ny nz]);
    f=reshape(S0*Ainv(2,:).',[nx ny nz]);
end