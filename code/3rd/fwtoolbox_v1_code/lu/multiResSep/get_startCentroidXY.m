function [cx, cy] = get_startCentroidXY(mask)

% find the center of the mass in the mask
im_mag = mask;
tx = sum(im_mag, 1).*(1:size(im_mag, 2));
ty = sum(im_mag, 2).*(1:size(im_mag, 1))';
cx = round(sum(tx)/sum(im_mag(:)));
cy = round(sum(ty)/sum(im_mag(:)));

