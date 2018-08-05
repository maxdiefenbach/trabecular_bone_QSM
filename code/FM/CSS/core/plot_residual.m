clear; close all; clc;

addpath(genpath('~/programs/BMRR/Postprocessing/CSS/'))

nTE = 6
dTE_s = 1e-3
TEmin_s = 1e-3
TE_s = TEmin_s + 0:(nTE-1) * dTE_s


fatfraction_percent = 30
deshielding_ppm = [0.9, 1.3, 1.6, 2.02, 2.24, 2.75, 4.2] - 4.7
relamps_percent = [8.9, 59.6, 5.9, 8.0, 5.9, 1.4, 3.9]
relamps_percent = 100 * relamps_percent / sum(relamps_percent)
fieldmap_Hz = 100
R2s_Hz = 0


sigamp = 100;
phases = 0;
centerfreq_Hz = 128e6;

pm = ones(length(deshielding_ppm)+1, 4)
pm(1, 1) = sigamp * (1 - fatfraction_percent/100)
pm(2:end, 1) = sigamp * fatfraction_percent/100 * relamps_percent/100;
pm(:, 2) = phases;
pm(1, 3) = fieldmap_Hz;
pm(2:end, 3) = fieldmap_Hz + centerfreq_Hz * deshielding_ppm * 1e-6;
pm(:, 4) = R2s_Hz;
pm

A = get_CSS_Amatrix(TE_s, pm)

s = build_signal(TE_s, pm)
