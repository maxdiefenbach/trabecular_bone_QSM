close all; clear all; clc;


I = ImDataParamsBMRR('/Users/maxdiefenbach/Projects/ImprovedWFI/data/reconstructed/Ankle/20161104/20161104_161845_1302_ImDataParams.mat')
magnetInhom_T = I.B0params.magnetInhom_T;
I.load_B0params

nii = make_nii(I.B0params.magnetInhom_T, ...
               I.ImDataParams.voxelSize_mm, ...
               I.ImDataParams.isoCenter_REC, ...
               [], ...
               I.fileID)
save_nii(nii, ['./' I.fileID '_magnetInhom.nii'])

I.save_array2nii(magnetInhom_T)