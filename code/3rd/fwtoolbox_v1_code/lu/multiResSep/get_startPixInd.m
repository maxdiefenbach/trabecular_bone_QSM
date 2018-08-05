function [start_pixel_ind, start_pixel_psi] = get_startPixInd(pd, mask)
% step 1
% look for the median 15 super-pixels
pdnzInd = find(mask~=0);
pdnz = pd(pdnzInd);
dum = abs(pdnz-median(pdnz)*ones(size(pdnz)));
[dum, sInd] = sort(dum);
% the indices of the 15 super-pixels
median_pixel_ind = pdnzInd(sInd(1:min(15, length(sInd))));

% step 2
% find the center of the mass in the mask
im_mag = mask;
tx = sum(im_mag, 1).*(1:size(im_mag, 2));
ty = sum(im_mag, 2).*(1:size(im_mag, 1))';
cx = round(sum(tx)/sum(im_mag(:)));
cy = round(sum(ty)/sum(im_mag(:)));

% step 3
% locate the super-pixel in the median group which is the closest to the
% center for the mass
[ht, wd] = size(pd);
my = mod(median_pixel_ind, ht); if my==0, my=ht; end
mx = floor(median_pixel_ind/ht);if mx==0, mx=1; end
dist2center = (mx-cx).^2+(my-cy).^2;
[dum, sInd] = min(dist2center);
start_pixel_ind = median_pixel_ind(sInd);
start_pixel_psi = pd(start_pixel_ind);