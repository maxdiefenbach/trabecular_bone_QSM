clear; close all; clc;

datadir = '../../../../data/reconstructed/'
subjectID = 'subject02'
iSlice = 88

%% load
subjectDir = fullfile(datadir, subjectID);

tmp = regexpdir(subjectDir, '\d{8}_\d{6}_\d{4}_echoMIP.nii.gz');
fmip = tmp{1}
nii_echoMIP = load_untouch_nii(fmip);
echoMIP = nii_echoMIP.img;

[~, fname] = fileparts(fmip);
fileID = fname(1:20)

fchi_close = fullfile(subjectDir, [fileID '_chi_closedFormL2_LBV_unwrap_CSS.nii.gz'])
nii_chi_close = load_untouch_nii(fchi_close);
chi_close = nii_chi_close.img;

fchi_cg = fullfile(subjectDir, [fileID '_chi_MEDIl2TVcg_LBV_unwrap_CSS.nii.gz'])
nii_chi_cg = load_untouch_nii(fchi_cg);
chi_cg = nii_chi_cg.img;

fchi_nesta = fullfile(subjectDir, [fileID '_chi_MEDIl1TVnesta_LBV_unwrap_CSS.nii.gz'])
nii_chi_nesta = load_untouch_nii(fchi_nesta);
chi_nesta = nii_chi_nesta.img;

fr2s = fullfile(subjectDir, [fileID '_R2s_CSS.nii.gz'])
nii_r2s = load_untouch_nii(fr2s);
r2s = nii_r2s.img;

fpdff = fullfile(subjectDir, [fileID '_PDFF_CSS.nii.gz'])
nii_pdff = load_untouch_nii(fpdff);
pdff = nii_pdff.img;


mask = echoMIP > 0.03 * max(echoMIP(:));
pdff = mask .* pdff;
r2s = mask .* r2s;

%% check
close all;
imagine(echoMIP, ...
        pdff, 'Window', [0, 100], ...
        r2s, 'Window', [0, 300], ...
        chi_nesta, 'Window', [-2, 2])

ix = 230
iy = 230

h = plot_orientations(echoMIP, ix, iy, iSlice)
set(h, 'Colormap', gray)
matlab2tikz('mip.tex', 'StandAlone', true)


pdff_range = [0, 100];
h = plot_orientations(pdff, ix, iy, iSlice)
axes = get(h, 'Children');
set(h, 'Colormap', inferno)
set(axes(1), 'CLim', pdff_range);
set(axes(2), 'CLim', pdff_range);
set(axes(3), 'CLim', pdff_range);
colorbar
matlab2tikz('pdff.tex', 'StandAlone', true)


r2s_range = [0, 300];
h = plot_orientations(r2s, ix, iy, iSlice)
axes = get(h, 'Children');
set(h, 'Colormap', magma)
set(axes(1), 'CLim', r2s_range);
set(axes(2), 'CLim', r2s_range);
set(axes(3), 'CLim', r2s_range);
colorbar
matlab2tikz('r2s.tex', 'StandAlone', true)


chi_range = [-2, 2];
h = plot_orientations(chi_nesta, ix, iy, iSlice)
axes = get(h, 'Children');
set(h, 'Colormap', viridis)
set(axes(1), 'CLim', chi_range);
set(axes(2), 'CLim', chi_range);
set(axes(3), 'CLim', chi_range);
colorbar
matlab2tikz('chi.tex', 'StandAlone', true)


close all;
figure()
colormap('viridis')
imagesc(chi_nesta(:, :, iSlice))
[Cx, Cy, chi_improfile] = improfile
r2s_improfile = improfile(r2s(:, :, iSlice), Cx, Cy)
pdff_improfile = improfile(pdff(:, :, iSlice), Cx, Cy)

save('improfiles.mat', 'Cx', 'Cy', 'chi_improfile', 'r2s_improfile', 'pdff_improfile')

T = table(Cx, Cy, pdff_improfile, r2s_improfile, chi_improfile)
writetable(T, 'improfiles.csv')

figure
subplot(1, 3, 1)
plot(pdff_improfile)
subplot(1, 3, 2)
plot(r2s_improfile)
subplot(1, 3, 3)
plot(chi_improfile)
matlab2tikz('improfiles.tex', 'StandAlone', true)