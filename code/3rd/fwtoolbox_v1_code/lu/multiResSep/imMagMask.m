function [mask, im1, im2, im3] = imMagMask(im1, im2, im3, magLevel)
if nargin==1,
    imall = im1;
    im1 = imall(:,:,1); im2 = imall(:,:,2); im3 = imall(:,:,3);
    imall = [im1(:)'; im2(:)'; im3(:)'];
end
if nargin==2,
    magLevel = im2;
if size(im1, 3) == 3,
    imall = im1;
else
    imall = repmat(im1(:,:,1), [1 1 3]);
end
    im1 = imall(:,:,1); im2 = imall(:,:,2); im3 = imall(:,:,3);
    imall = [im1(:)'; im2(:)'; im3(:)'];
end
if nargin>=3,
    imall = [im1(:)'; im2(:)'; im3(:)'];
end
% define a mask, which blocks the background (i.e., signal
% magnitude<magLevel*range)
if ~exist('magLevel')
magLevel = .1/2;
end
imall_mag = sum(abs(imall), 1);
mask = zeros(size(im1));
mask(find(imall_mag>magLevel.*range(imall_mag))) = 1;
mask = im2bw(mask);
% figure; imshow(mask); drawnow;

if nargout > 1
    im1 = im1.*mask;
    im2 = im2.*mask;
    im3 = im3.*mask;
end
