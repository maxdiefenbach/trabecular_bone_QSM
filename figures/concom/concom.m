clear; close all; clc;

datadir = '../../../../data/reconstructed/'
subjectID = 'subject04'
iSlice = 80

%% load
subjectDir = fullfile(datadir, subjectID);

tmp = regexpdir(subjectDir, '\d{8}_\d{6}_\d{4}_echoMIP.nii.gz');
fmip = tmp{1}
nii_echoMIP = load_untouch_nii(fmip);
echoMIP = nii_echoMIP.img;

[~, fname] = fileparts(fmip);
fileID = fname(1:20)

% WFI
fpdff = fullfile(subjectDir, [fileID '_PDFF_CSS.nii.gz'])
nii_pdff = load_untouch_nii(fpdff);
pdff = nii_pdff.img;

fr2s = fullfile(subjectDir, [fileID '_R2s_CSS.nii.gz'])
nii_r2s = load_untouch_nii(fr2s);
r2s = nii_r2s.img;

% LBV
flbv = fullfile(subjectDir, [fileID '_localRDF_LBV_unwrap_CSS.nii.gz'])
nii_lbv = load_untouch_nii(flbv);
lbv = nii_lbv.img;

flbv_close = fullfile(subjectDir, [fileID '_QSMparams_closedFormL2_CSSLBV_unwrap.nii.gz'])
nii_lbv_close = load_untouch_nii(flbv_close);
lbv_close = nii_lbv_close.img;

flbv_cg = fullfile(subjectDir, [fileID '_QSMparams_MEDIl2TVcg_CSSLBV_unwrap.nii.gz'])
nii_lbv_cg = load_untouch_nii(flbv_cg);
lbv_cg = nii_lbv_cg.img;

flbv_nesta = fullfile(subjectDir, [fileID '_QSMparams_MEDIl1TVnesta_CSSLBV_unwrap.nii.gz'])
nii_lbv_nesta = load_untouch_nii(flbv_nesta);
lbv_nesta = nii_lbv_nesta.img;

% PDF
fpdf = fullfile(subjectDir, [fileID '_localRDF_PDF_unwrap_CSS.nii.gz'])
nii_pdf = load_untouch_nii(fpdf);
pdf = nii_pdf.img;

fpdf_close = fullfile(subjectDir, [fileID '_QSMparams_closedFormL2_CSSPDF_unwrap.nii.gz'])
nii_pdf_close = load_untouch_nii(fpdf_close);
pdf_close = nii_pdf_close.img;

fpdf_cg = fullfile(subjectDir, [fileID '_QSMparams_MEDIl2TVcg_CSSPDF_unwrap.nii.gz'])
nii_pdf_cg = load_untouch_nii(fpdf_cg);
pdf_cg = nii_pdf_cg.img;

fpdf_nesta = fullfile(subjectDir, [fileID '_QSMparams_MEDIl1TVnesta_CSSPDF_unwrap.nii.gz'])
nii_pdf_nesta = load_untouch_nii(fpdf_nesta);
pdf_nesta = nii_pdf_nesta.img;

mask = echoMIP >= 0.1 * max(echoMIP(:));
mask = fill_mask3D(mask, 1);
pdf = mask .* pdf;


% plots
close all;
rdf_range = [-1, 1]
chi_range = [-2, 2]
imagine(lbv, 'Window', rdf_range, ...
        lbv_close, 'Window', chi_range, ...
        lbv_cg, 'Window', chi_range, ...
        lbv_nesta, 'Window', chi_range)

imagine(pdf, 'Window', rdf_range, ...
        pdf_close, 'Window', chi_range, ...
        pdf_cg, 'Window', chi_range, ...
        pdf_nesta, 'Window', chi_range)


close all;
ix = 230
iy = 230
iz = iSlice

h = plot_orientations(echoMIP, ix, iy, iz)
set(h, 'Colormap', gray)
matlab2tikz('echoMIP.tex', 'StandAlone', true)

% lbv
h = plot_orientations(lbv, ix, iy, iz)
set(h, 'Colormap', plasma)
axes = get(h, 'Children')
set(axes(1), 'CLim', rdf_range)
set(axes(2), 'CLim', rdf_range)
set(axes(3), 'CLim', rdf_range)
colorbar('peer', axes(1))
matlab2tikz('lbv.tex', 'StandAlone', true)

h = plot_orientations(lbv_close, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('lbv_close.tex', 'StandAlone', true)

h = plot_orientations(lbv_cg, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('lbv_cg.tex', 'StandAlone', true)

h = plot_orientations(lbv_nesta, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('lbv_nesta.tex', 'StandAlone', true)

% pdf
h = plot_orientations(pdf, ix, iy, iz)
set(h, 'Colormap', plasma)
axes = get(h, 'Children')
set(axes(1), 'CLim', rdf_range)
set(axes(2), 'CLim', rdf_range)
set(axes(3), 'CLim', rdf_range)
colorbar('peer', axes(1))
matlab2tikz('pdf.tex', 'StandAlone', true)

h = plot_orientations(pdf_close, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('pdf_close.tex', 'StandAlone', true)

h = plot_orientations(pdf_cg, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('pdf_cg.tex', 'StandAlone', true)

h = plot_orientations(pdf_nesta, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('pdf_nesta.tex', 'StandAlone', true)


%% concom

% WFI
fpdff = fullfile(subjectDir, [fileID '_PDFF_CSS_concomCorr.nii.gz'])
nii_pdff = load_untouch_nii(fpdff);
pdff = nii_pdff.img;

fr2s = fullfile(subjectDir, [fileID '_R2s_CSS_concomCorr.nii.gz'])
nii_r2s = load_untouch_nii(fr2s);
r2s = nii_r2s.img;

% LBV
flbv = fullfile(subjectDir, [fileID '_localRDF_LBV_unwrap_CSS_concomCorr.nii.gz'])
nii_lbv = load_untouch_nii(flbv);
lbv = nii_lbv.img;

flbv_close = fullfile(subjectDir, [fileID '_QSMparams_closedFormL2_CSS_concomCorrLBV_unwrap.nii.gz'])
nii_lbv_close = load_untouch_nii(flbv_close);
lbv_close = nii_lbv_close.img;

flbv_cg = fullfile(subjectDir, [fileID '_QSMparams_MEDIl2TVcg_CSS_concomCorrLBV_unwrap.nii.gz'])
nii_lbv_cg = load_untouch_nii(flbv_cg);
lbv_cg = nii_lbv_cg.img;

flbv_nesta = fullfile(subjectDir, [fileID '_QSMparams_MEDIl1TVnesta_CSS_concomCorrLBV_unwrap.nii.gz'])
nii_lbv_nesta = load_untouch_nii(flbv_nesta);
lbv_nesta = nii_lbv_nesta.img;

% PDF
fpdf = fullfile(subjectDir, [fileID '_localRDF_PDF_unwrap_CSS_concomCorr.nii.gz'])
nii_pdf = load_untouch_nii(fpdf);
pdf = nii_pdf.img;

fpdf_close = fullfile(subjectDir, [fileID '_QSMparams_closedFormL2_CSS_concomCorrPDF_unwrap.nii.gz'])
nii_pdf_close = load_untouch_nii(fpdf_close);
pdf_close = nii_pdf_close.img;

fpdf_cg = fullfile(subjectDir, [fileID '_QSMparams_MEDIl2TVcg_CSS_concomCorrPDF_unwrap.nii.gz'])
nii_pdf_cg = load_untouch_nii(fpdf_cg);
pdf_cg = nii_pdf_cg.img;

fpdf_nesta = fullfile(subjectDir, [fileID '_QSMparams_MEDIl1TVnesta_CSS_concomCorrPDF_unwrap.nii.gz'])
nii_pdf_nesta = load_untouch_nii(fpdf_nesta);
pdf_nesta = nii_pdf_nesta.img;

mask = echoMIP >= 0.1 * max(echoMIP(:));
mask = fill_mask3D(mask, 1);
pdf = mask .* pdf;


close all;
ix = 230
iy = 230
iz = iSlice

% lbv
h = plot_orientations(lbv, ix, iy, iz)
set(h, 'Colormap', plasma)
axes = get(h, 'Children')
set(axes(1), 'CLim', rdf_range)
set(axes(2), 'CLim', rdf_range)
set(axes(3), 'CLim', rdf_range)
colorbar('peer', axes(1))
matlab2tikz('lbv_concom.tex', 'StandAlone', true)

h = plot_orientations(lbv_close, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('lbv_close_concom.tex', 'StandAlone', true)

h = plot_orientations(lbv_cg, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('lbv_cg_concom.tex', 'StandAlone', true)

h = plot_orientations(lbv_nesta, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('lbv_nesta_concom.tex', 'StandAlone', true)

% pdf
h = plot_orientations(pdf, ix, iy, iz)
set(h, 'Colormap', plasma)
axes = get(h, 'Children')
set(axes(1), 'CLim', rdf_range)
set(axes(2), 'CLim', rdf_range)
set(axes(3), 'CLim', rdf_range)
colorbar('peer', axes(1))
matlab2tikz('pdf_concom.tex', 'StandAlone', true)

h = plot_orientations(pdf_close, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('pdf_close_concom.tex', 'StandAlone', true)

h = plot_orientations(pdf_cg, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('pdf_cg_concom.tex', 'StandAlone', true)

h = plot_orientations(pdf_nesta, ix, iy, iz)
set(h, 'Colormap', viridis)
axes = get(h, 'Children')
set(axes(1), 'CLim', chi_range)
set(axes(2), 'CLim', chi_range)
set(axes(3), 'CLim', chi_range)
colorbar('peer', axes(1))
matlab2tikz('pdf_nesta_concom.tex', 'StandAlone', true)
