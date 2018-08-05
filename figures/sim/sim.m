clear; close all; clc;

datadir = '../../../../data/simulation/results/';

fileID = 'P4_B0dir001_nonUTE'


% chi
ptype = '_mask'
filename = [fileID ptype '.nii.gz'];
nii = load_nii(fullfile(datadir, filename))
center = ceil((size(nii.img) + 1)/2)
nii.img = padarray(nii.img, center);
iSlice = ceil((size(nii.img, 3) + 1)/2)
figure
colormap('gray')
imagesc(-nii.img(:, :, iSlice))
axis image
colorbar
matlab2tikz(strrep(filename, '.nii.gz', '.tex'), 'standalone', true)

% RDF
ptype = '_mask'
filename = [fileID ptype '.nii.gz'];
nii = load_nii(fullfile(datadir, filename))
center = ceil((size(nii.img) + 1)/2)
padsize = [center(1), center(2), 0]
nii.img = padarray(nii.img, padsize);
DataParams.chimap_ppm = nii.img(:, :, (center(3)-5):(center(3)+4));
DataParams.voxelSize_mm = [1, 1, 1];
DataParams.B0dir = [1, 0, 0];
RDF_ppm = forwardSimulate_RDF_ppm(DataParams);
iSlice = ceil((size(RDF_ppm, 3) + 1)/2)
figure
colormap('plasma')
imagesc(RDF_ppm(:, :, iSlice))
axis image
colorbar
matlab2tikz(strrep(strrep(filename, 'mask', 'RDF'), '.nii.gz', '.tex'), 'standalone', true)


% signal
figure
colormap('gray')
subplot(2, 2, 1)
imagesc(-DataParams.chimap_ppm(:, :, ceil((size(DataParams.chimap_ppm, 3)+1) / 2)))
axis image
subplot(2, 2, 2)
imagesc(-DataParams.chimap_ppm(:, :, ceil((size(DataParams.chimap_ppm, 3)+1) / 2)))
axis image
subplot(2, 2, 3)
imagesc(RDF_ppm(:, :, ceil((size(RDF_ppm, 3)+1) / 2)))
axis image
subplot(2, 2, 4)
imagesc(RDF_ppm(:, :, ceil((size(RDF_ppm, 3)+1) / 2)) * 100)
axis image
matlab2tikz(strrep(strrep(filename, 'mask', 'signal'), '.nii.gz', '.tex'), 'standalone', true)


ptype = '_R2s'
filename = [fileID ptype '.nii.gz'];
nii = load_nii(fullfile(datadir, filename))
iSlice = ceil((size(nii.img, 3) + 1)/2)
figure
colormap('magma')
imagesc(nii.img(:, :, iSlice))
axis image
caxis([0, 150])
colorbar
matlab2tikz(strrep(filename, '.nii.gz', '.tex'), 'standalone', true)

ptype = '_chi'
filename = [fileID ptype '.nii.gz'];
nii = load_nii(fullfile(datadir, filename))
iSlice = ceil((size(nii.img, 3) + 1)/2)
figure
colormap('viridis')
imagesc(nii.img(:, :, iSlice))
axis image
colorbar
matlab2tikz(strrep(filename, '.nii.gz', '.tex'), 'standalone', true)

ptype = '_RDF_downsamp'
filename = [fileID ptype '.nii.gz'];
nii = load_nii(fullfile(datadir, filename))
iSlice = ceil((size(nii.img, 3) + 1)/2)
figure
colormap('plasma')
imagesc(nii.img(:, :, iSlice))
axis image
colorbar
matlab2tikz(strrep(filename, '.nii.gz', '.tex'), 'standalone', true)


ptype = '_signal_downsamp'
filename = [fileID ptype '.nii.gz'];
nii = load_untouch_nii(fullfile(datadir, filename))
% nii.img = permute(nii.img, [2, 1, 3]);
iSlice = ceil((size(nii.img, 3) + 1)/2)
figure
colormap('gray')
subplot(2, 2, 1)
imagesc(abs(nii.img(:, :, iSlice, 1)))
axis image
subplot(2, 2, 2)
imagesc(angle(nii.img(:, :, iSlice, 2)))
axis image
subplot(2, 2, 3)
imagesc(abs(nii.img(:, :, iSlice, 1)))
axis image
subplot(2, 2, 4)
imagesc(angle(nii.img(:, :, iSlice, 2)))
axis image
matlab2tikz(strrep(filename, '.nii.gz', '.tex'), 'standalone', true)
