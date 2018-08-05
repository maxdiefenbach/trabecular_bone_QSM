clear; close all; clc;

datadir = '../../../../data/reconstructed/'
subjectID = 'subject01'
iSlice = 66

%% load
subjectDir = fullfile(datadir, subjectID);

tmp = regexpdir(subjectDir, '\d{8}_\d{6}_\d{4}_echoMIP.nii.gz');
fmip = tmp{1}
nii_echoMIP = load_untouch_nii(fmip);
echoMIP = nii_echoMIP.img;

[~, fname] = fileparts(fmip);
fileID = fname(1:20)

fchi_close = fullfile(subjectDir, [fileID '_QSMparams_closedFormL2_CSSLBV_unwrap.nii.gz'])
nii_chi_close = load_untouch_nii(fchi_close);
chi_close = nii_chi_close.img;

fchi_cg = fullfile(subjectDir, [fileID '_QSMparams_MEDIl2TVcg_CSSLBV_unwrap.nii.gz'])
nii_chi_cg = load_untouch_nii(fchi_cg);
chi_cg = nii_chi_cg.img;

fchi_nesta = fullfile(subjectDir, [fileID '_QSMparams_MEDIl1TVnesta_CSSLBV_unwrap.nii.gz'])
nii_chi_nesta = load_untouch_nii(fchi_nesta);
chi_nesta = nii_chi_nesta.img;

fr2s = fullfile(subjectDir, [fileID '_R2s_CSS.nii.gz'])
nii_r2s = load_untouch_nii(fr2s);
r2s = nii_r2s.img;

fpdff = fullfile(subjectDir, [fileID '_PDFF_CSS.nii.gz'])
nii_pdff = load_untouch_nii(fpdff);
pdff = nii_pdff.img;

fbssfp = fullfile(subjectDir, [fileID '_bssfp_N4biasCorr_resamp.nii.gz'])
nii_bssfp = load_untouch_nii(fbssfp);
bssfp = nii_bssfp.img;
bssfp = flip(flipud(permute(nii_bssfp.img, [2, 1, 3])), 3);

%% check
chi_range = [-2, 2]
close all;
imagine(chi_close, 'Window', chi_range, ...
        chi_cg, 'Window', chi_range, ...
        chi_nesta, 'Window', chi_range,...
        'Panels', [1, 3])
imagine(echoMIP)
imagine(bssfp)


%% plot

% anatomical images
figure
colormap('gray')
subplot(1, 2, 1)
imagesc(echoMIP(:, :, iSlice))
axis image
axis off

subplot(1, 2, 2)
imagesc(bssfp(:, :, 35))
axis image
axis off
matlab2tikz('anatomical_images.tex', 'StandAlone', true)

% pdff
figure
colormap('inferno')
imagesc(pdff(:, :, iSlice))
caxis([0, 100])
axis image
axis off
colorbar
matlab2tikz('pdff.tex', 'StandAlone', true)

% r2s
figure
colormap('magma')
imagesc(pdff(:, :, iSlice))
caxis([0, 100])
axis image
axis off
colorbar
matlab2tikz('r2s.tex', 'StandAlone', true)

% chi images
figure
colormap('viridis')
nrow = 1;
ncol = 3;

subplot(nrow, ncol, 1)
imagesc(chi_close(:, :, iSlice))
caxis(chi_range)
axis image
axis off

subplot(nrow, ncol, 2)
imagesc(chi_cg(:, :, iSlice))
caxis(chi_range)
axis image
axis off

subplot(nrow, ncol, 3)
imagesc(chi_nesta(:, :, iSlice))
caxis(chi_range)
axis image
axis off
colorbar

matlab2tikz('chi_images.tex', 'StandAlone', true)
