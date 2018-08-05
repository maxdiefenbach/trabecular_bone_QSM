clear; close all; clc;

datadir = '../../../../data/reconstructed/'
subjectID = 'subject05'
iSlice = 68

%% load
subjectDir = fullfile(datadir, subjectID);

% src
tmp = regexpdir(subjectDir, '\d{8}_\d{6}_\d{4}_ImDataParams.mat');
filename = tmp{1}
I = ImDataParamsBMRR(filename)
I.plot

mask3D = I.get_tissueMask(5);
mask4D = repmat(mask3D, [1, 1, 1, 9]);

P = mask4D .* I.phase;
M = mask4D .* I.magnitude;


close all;

figure
colormap('gray')

subplot(3, 2, 1)
imagesc(M(:, :, iSlice, 1))
axis image
subplot(3, 2, 2)
imagesc(P(:, :, iSlice, 1))
axis image

subplot(3, 2, 3)
imagesc(M(:, :, iSlice, 2))
axis image
subplot(3, 2, 4)
imagesc(P(:, :, iSlice, 2))
axis image

subplot(3, 2, 5)
imagesc(M(:, :, iSlice, end))
axis image
subplot(3, 2, 6)
imagesc(P(:, :, iSlice, end))
axis image

matlab2tikz('src.tex', 'standalone', true)


% WFI
tmp = regexpdir(subjectDir, '\d{8}_\d{6}_\d{4}_WFIparams_graphcut.mat');
filename = tmp{1}
I.load_WFIparams(filename)

I.WFIparams
I.set_fatFraction_percent

close all;
figure
colormap('inferno')
imagesc(mask3D(:, :, iSlice) .* I.WFIparams.fatFraction_percent(:, :, iSlice))
caxis([0, 100])
axis image
colorbar
matlab2tikz('pdff.tex', 'standalone', true)

figure
colormap('magma')
imagesc(mask3D(:, :, iSlice) .* I.WFIparams.R2s_Hz(:, :, iSlice))
caxis([0, 300])
colorbar
axis image
matlab2tikz('R2s.tex', 'standalone', true)

figure
colormap('plasma')
imagesc(mask3D(:, :, iSlice) .* I.WFIparams.fieldmap_Hz(:, :, iSlice))
caxis([-300, 300])
colorbar
axis image
matlab2tikz('total_field_Hz.tex', 'standalone', true)


% BFR
tmp = regexpdir(subjectDir, '\d{8}_\d{6}_\d{4}_BFRparamsLBV_graphcut.mat');
filename = tmp{1}
I.load_BFRparams(filename)
I.BFRparams


figure
colormap('plasma')
imagesc(I.BFRparams.totalRDF_ppm(:, :, iSlice))
colorbar
axis image
matlab2tikz('total_field_ppm.tex', 'standalone', true)


figure
colormap('plasma')
imagesc(I.BFRparams.localRDF_ppm(:, :, iSlice))
colorbar
axis image
matlab2tikz('local_field_ppm.tex', 'standalone', true)


% QSM
tmp = regexpdir(subjectDir, '\d{8}_\d{6}_\d{4}_QSMparams_MEDIl1TVnesta_CSSLBV_unwrap.mat');
filename = tmp{1}
I.load_QSMparams(filename)


figure
colormap('viridis')
imagesc(mask3D(:, :, iSlice) .* I.QSMparams.chimap_ppm(:, :, iSlice))
caxis([-2, 2])
colorbar
axis image
matlab2tikz('chi.tex', 'standalone', true)
