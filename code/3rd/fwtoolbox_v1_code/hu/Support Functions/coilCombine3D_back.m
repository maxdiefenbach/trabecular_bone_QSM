% Function: coilCombine3D
%
% Description: combine multi-coil image sequences
%
% Based on: Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of
% phased array MR imagery. Magn Reson Med 2000;43:682-690
% 
% Parameters:
% im1: the multi-coil images (size [nx,ny,nz,ncoils,nTE])
%
% Returns: 
% im2: the coil-combined images (size [nx,ny,nz,nTE])
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: November 21, 2011

function im2 = coilCombine3D(im1)
% sz = size(im1,3);
[sx,sy,sz,~,sTE]=size(im1);
im2 = zeros(sx,sy,sz,1,sTE);

for kz=1:sz
  im2(:,:,kz,1,:) = coilCombine(im1(:,:,kz,:,:));
end
