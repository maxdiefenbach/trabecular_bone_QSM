clear; close all; clc;

% setup
p = path;                               % for later tear down
restoredefaultpath;
addpath(genpath('../'))



filename = '../data/20180406_144018_0302_ImDataParams.mat'
I = ImDataParamsBMRR(filename)

I.load_WFIparams('../data/20180406_144018_0302_WFIparams_CSS_unwrap.mat')
% I = WFI_subject(filename)
I.plot_WFIparams

I.load_BFRparams('../data/20180406_144018_0302_BFRparams_LBV_CSS_unwrap.mat')
% I = BFR_subject(filename)
I.plot_BFRparams

I.load_QSMparams('../data/20180406_144018_0302_QSMparams_MEDIl1TVnesta_LBV_CSS_unwrap.mat')
% I = QSM_subject(filename)
I.plot_QSMparams

echoMIP = I.get_echoMIP;
mask = echoMIP > (0.05 * max(echoMIP(:)));
iz = ceil((size(echoMIP, 3) + 1) / 2)


close all;
h = figure
ax = subplot(2, 3, 1)
colormap(ax, gray)
imagesc(mask(:, :, iz) .* echoMIP(:, :, iz))
xlabel('MIP_{TE} [a.u.]')
set(gca,'xtick',[])
set(gca,'ytick',[])

ax = subplot(2, 3, 2)
imagesc(mask(:, :, iz) .* I.WFIparams.fatFraction_percent(:, :, iz))
colormap(ax, inferno)
caxis([0, 100])
colorbar
xlabel('PDFF [%]')
set(gca,'xtick',[])
set(gca,'ytick',[])

ax = subplot(2, 3, 3)
imagesc(mask(:, :, iz) .* I.WFIparams.R2s_Hz(:, :, iz))
colormap(ax, magma)
caxis([0, 300])
colorbar
xlabel('R_2^*')
set(gca,'xtick',[])
set(gca,'ytick',[])

ax = subplot(2, 3, 4)
imagesc(mask(:, :, iz) .* I.WFIparams.fieldmap_Hz_unwrap(:, :, iz))
colormap(ax, plasma)
colorbar
xlabel('field map [Hz]')
set(gca,'xtick',[])
set(gca,'ytick',[])

ax = subplot(2, 3, 5)
colormap(ax, plasma)
imagesc(mask(:, :, iz) .* I.BFRparams.localRDF_ppm(:, :, iz))
caxis([-1, 1])
colorbar
xlabel('local field [ppm]')
set(gca,'xtick',[])
set(gca,'ytick',[])

ax = subplot(2, 3, 6)
imagesc(mask(:, :, iz) .* I.QSMparams.chimap_ppm(:, :, iz))
colormap(ax, viridis)
caxis([-2, 2])
colorbar
xlabel('susceptibility [ppm]')
set(gca,'xtick',[])
set(gca,'ytick',[])

saveas(h, 'output.png')


% teardown: restore previous matlab path
path(p)