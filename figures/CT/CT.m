clear; close all; clc;

datadir = '../../../../data/reconstructed/'

subjectID = 'subjectCT16'

subjectDir = fullfile(datadir, subjectID);

%% load
tmp = regexpdir(subjectDir, '\d{8}_\d{6}_\d{4}_echoMIP.nii.gz');
fmip = tmp{1}
nii_mip = load_untouch_nii(fmip);
mip = nii_mip.img;

fchi = strrep(fmip, 'echoMIP', 'chi_MEDIl1TVnesta_LBV_CSS')
nii_chi = load_untouch_nii(fchi);
chi = nii_chi.img;

fr2s = strrep(fmip, 'echoMIP', 'R2s_CSS')
nii_r2s = load_untouch_nii(fr2s);
r2s = nii_r2s.img;

tmp = regexpdir(fullfile(subjectDir, 'CT'), '.*_rightFoot_resamp2echoMIP.nii.gz');
fct = tmp{1}
nii_ct = load_untouch_nii(fct);
ct = nii_ct.img;

%% plot
ix = 277
iy = 296
iz = 61

close all;
h = plot_orientations(mip, ix, iy, iz)
set(h, 'Colormap', gray)
matlab2tikz('mip1.tex', 'StandAlone', true)

h = plot_orientations(ct, ix, iy, iz)
set(h, 'Colormap', gray)
colorbar
matlab2tikz('ct1.tex', 'StandAlone', true)

h = plot_orientations(chi, ix, iy, iz)
set(h, 'Colormap', viridis)
chi_range = [-2, 2];
axes = get(h, 'Children');
set(axes(1), 'CLim', chi_range);
set(axes(2), 'CLim', chi_range);
set(axes(3), 'CLim', chi_range);
colorbar
matlab2tikz('chi1.tex', 'StandAlone', true)

h = plot_orientations(r2s, ix, iy, iz)
set(h, 'Colormap', magma)
axes = get(h, 'Children');
r2s_range = [0, 300];
set(axes(1), 'CLim', r2s_range);
set(axes(2), 'CLim', r2s_range);
set(axes(3), 'CLim', r2s_range);
colorbar
matlab2tikz('r2s1.tex', 'StandAlone', true)



subjectID = 'subjectCT17'

subjectDir = fullfile(datadir, subjectID);

%% load
tmp = regexpdir(subjectDir, '\d{8}_\d{6}_\d{4}_echoMIP.nii.gz');
fmip = tmp{1}
nii_mip = load_untouch_nii(fmip);
mip = nii_mip.img;

fchi = strrep(fmip, 'echoMIP', 'chi_MEDIl1TVnesta_LBV_CSS')
nii_chi = load_untouch_nii(fchi);
chi = nii_chi.img;

fr2s = strrep(fmip, 'echoMIP', 'R2s_CSS')
nii_r2s = load_untouch_nii(fr2s);
r2s = nii_r2s.img;

tmp = regexpdir(fullfile(subjectDir, 'CT'), '.*_rightFoot_resamp2echoMIP.nii.gz');
fct = tmp{1}
nii_ct = load_untouch_nii(fct);
ct = nii_ct.img;

%% plot
ix = 247
iy = 282
iz = 81

close all;
h = plot_orientations(mip, ix, iy, iz)
set(h, 'Colormap', gray)
matlab2tikz('mip2.tex', 'StandAlone', true)

h = plot_orientations(ct, ix, iy, iz)
set(h, 'Colormap', gray)
colorbar
matlab2tikz('ct2.tex', 'StandAlone', true)

h = plot_orientations(chi, ix, iy, iz)
set(h, 'Colormap', viridis)
chi_range = [-1, 1];
axes = get(h, 'Children');
set(axes(1), 'CLim', chi_range);
set(axes(2), 'CLim', chi_range);
set(axes(3), 'CLim', chi_range);
colorbar
matlab2tikz('chi2.tex', 'StandAlone', true)

h = plot_orientations(r2s, ix, iy, iz)
set(h, 'Colormap', magma)
axes = get(h, 'Children');
r2s_range = [0, 300];
set(axes(1), 'CLim', r2s_range);
set(axes(2), 'CLim', r2s_range);
set(axes(3), 'CLim', r2s_range);
colorbar
matlab2tikz('r2s2.tex', 'StandAlone', true)