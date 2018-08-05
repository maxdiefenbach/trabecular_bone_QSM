%% Function name: getMagnitudeWeight
%%
%% Description: Calculates a "magnitude weight"-image according to
%%              Eqs. 5-6 in:
%%
%% Berglund J, Johansson L, Ahlström H, Kullberg J. Three-point Dixon method enables whole-body 
%% water and fat imaging of obese subjects. Magn Reson Med. 2010, Jun;63(6):1659-1668.
%%
%% Input:
%%   - S: complex multi-echo image, size[nx,ny,nz,N]
%%
%% Output:
%%   - mw: real-valued magnitude weight image, size[nx,ny,nz]
%%
%% Author: Johan Berglund
%% Date created: November 16, 2011
%% Date last modified: November 18, 2011

function mw = getMagnitudeWeight(S)
    %% Calculate magnitude according to Eq. 5:
    m=sum(abs(S),4);

    %% Find threshold using Otsu's method:
    m=m/max(reshape(m,[],1)); %Normalize magnitude
    mt = graythresh(m);                                                         

    %% Find mean value of foreground voxels:
    mf = mean(mean(m(m>mt)));

    %% Calculate magnitude weight according to Eq. 6:
    mw=1.0./(1+exp(3*(m-mt)/(mt-mf)));
end