clear; close all; clc;


subjectID = 'subject01'
WFImethod = 'CSS'
BFRmethod = 'PDF_unwrap'

datadir = '../../../../data/reconstructed/lcurve'


filename = [subjectID '_lcurve_close_' WFImethod '_' BFRmethod '.mat']
load(fullfile(datadir, filename))

filename = [subjectID '_lcurve_cg_' WFImethod '_' BFRmethod '.mat']
load(fullfile(datadir, filename))

filename = [subjectID '_lcurve_nesta_' WFImethod '_' BFRmethod '.mat']
load(fullfile(datadir, filename))

tmp = regexpdir(['../../../../data/reconstructed/' subjectID], '.*_echoMIP.nii.gz')
nii = load_untouch_nii(tmp{1})
echoMIP = nii.img;
mask = fill_mask3D(echoMIP >= 0.05 * max(echoMIP(:)), 1);
imagine(mask)


% close all;
% figure
% loglog(regparamList, resnormList_close, 'o-')
% figure
% loglog(regparamList, resnormList_cg, 'o-')
% figure
% loglog(regparamList, resnormList_nesta, 'o-')

close all;

imagine(echoMIP)

ix = 234
iy = 217
iz = 71

close all;
% closed form
regind = 4;
h = plot_orientations(mask .* squeeze(chimaps_close(:, :, :, regind)), ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
chi_range = [-2, 2];
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz([subjectID '_chi_close_' WFImethod '_' BFRmethod '_' num2str(regind) '.tex'], 'standalone', true)

regind = 6;
h = plot_orientations(mask .* squeeze(chimaps_close(:, :, :, regind)), ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
chi_range = [-2, 2];
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz([subjectID '_chi_close_' WFImethod '_' BFRmethod '_' num2str(regind) '.tex'], 'standalone', true)

regind = 11;
h = plot_orientations(mask .* squeeze(chimaps_close(:, :, :, regind)), ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
chi_range = [-2, 2];
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz([subjectID '_chi_close_' WFImethod '_' BFRmethod '_' num2str(regind) '.tex'], 'standalone', true)



% cg
regind = 4;
h = plot_orientations(mask .* squeeze(chimaps_cg(:, :, :, regind)), ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
chi_range = [-2, 2];
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz([subjectID '_chi_cg_' WFImethod '_' BFRmethod '_' num2str(regind) '.tex'], 'standalone', true)

regind = 6;
h = plot_orientations(mask .* squeeze(chimaps_cg(:, :, :, regind)), ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
chi_range = [-2, 2];
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz([subjectID '_chi_cg_' WFImethod '_' BFRmethod '_' num2str(regind) '.tex'], 'standalone', true)

regind = 11;
h = plot_orientations(mask .* squeeze(chimaps_cg(:, :, :, regind)), ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
chi_range = [-2, 2];
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz([subjectID '_chi_cg_' WFImethod '_' BFRmethod '_' num2str(regind) '.tex'], 'standalone', true)




% nesta
regind = 4;
h = plot_orientations(mask .* squeeze(chimaps_nesta(:, :, :, regind)), ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
chi_range = [-2, 2];
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz([subjectID '_chi_nesta_' WFImethod '_' BFRmethod '_' num2str(regind) '.tex'], 'standalone', true)

regind = 6;
h = plot_orientations(mask .* squeeze(chimaps_nesta(:, :, :, regind)), ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
chi_range = [-2, 2];
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz([subjectID '_chi_nesta_' WFImethod '_' BFRmethod '_' num2str(regind) '.tex'], 'standalone', true)

regind = 11;
h = plot_orientations(mask .* squeeze(chimaps_nesta(:, :, :, regind)), ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
chi_range = [-2, 2];
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz([subjectID '_chi_nesta_' WFImethod '_' BFRmethod '_' num2str(regind) '.tex'], 'standalone', true)


close all;
imagine(chimaps_close, 'Window', [-2, 2])
imagine(chimaps_cg, 'Window', [-2, 2])
imagine(chimaps_nesta, 'Window', [-2, 2])

figure
loglog(regparamList, resnormList_close, 'o-')
hold on
loglog(regparamList, resnormList_cg, 'o-')
loglog(regparamList, resnormList_nesta, 'o-')
ylim([1e3, 1e6])
ylabel('discrepancy')
xlabel('regularization parameter')
legend({'l2-TV cl form sol.', 'l2-TV conj. gradients', 'l1-TV NESTA'})
matlab2tikz([subjectID '_lcurve_plot.tex'], 'standalone', true)
