% 
% example: water and fat complex image separation
%

clear; close all;

echo on;
% First example
echo off;

load exampleData;
wat = abs(wat);
fat = abs(fat);

df_wf = -220; 
TE = TE./1000; % Convert to seconds.

% synthesize the data
im1 = wat + fat.*exp(i*2*pi*df_wf*TE(1)); 
im2 = wat + fat.*exp(i*2*pi*df_wf*TE(2)); 
im3 = wat + fat.*exp(i*2*pi*df_wf*TE(3)); 

% fmap = fmap + randn(size(fmap)).*10;
im1 = im1.*exp(i*2*pi*fmap.*TE(1));
im2 = im2.*exp(i*2*pi*fmap.*TE(2));
im3 = im3.*exp(i*2*pi*fmap.*TE(3));
im = cat(3, im1, im2, im3);

% psi = est_fieldMap(TE, im, 0);
% [w, f, sc, resd] = lsSep_wf(TE, im, psi);
% [w, f, psi] = idealScott(TE.*1000, im1, im2, im3);
[w, f, psi] = multiResSeparation(TE.*1000, im1, im2, im3);
pause(1);
echo on;
% Second example
echo off;
clear; 
load sample-data6;
[w, f, psi] = multiResSeparation(TE, im1, im2, im3);