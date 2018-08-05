import os
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt


def resample_bssfp(subjectID):
    datadir = '../../../../data/reconstructed'
    subjectDir = os.path.join(datadir, subjectID)
    fbssfp = [os.path.join(subjectDir, f) for f in os.listdir(subjectDir)
              if f.endswith('_bssfp_N4biasCorr.nii.gz')][0]
    print(subjectDir, os.listdir(subjectDir), fbssfp)
    fileID = os.path.basename(fbssfp)[0:20]
    ftrafo = os.path.join(subjectDir, fileID + '_transform_mip2bssfp.tfm')
    fmip = os.path.join(subjectDir, fileID + '_echoMIP.nii.gz')

    Image_bssfp = sitk.ReadImage(fbssfp)
    Image_mip = sitk.ReadImage(fmip)
    transform = sitk.ReadTransform(ftrafo).GetInverse()
    
    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(Image_mip)
    resample.SetTransform(transform)
    resample.SetInterpolator(sitk.sitkLinear)
    Image_bssfp_resamp = resample.Execute(Image_bssfp)

    fout = fbssfp.replace('_bssfp_N4biasCorr.nii.gz', '_bssfp_N4biasCorr_resamp.nii.gz')
    sitk.WriteImage(Image_bssfp_resamp, fout)
    print('Wrote {}.'.format(fout))

subjectID = 'subject01'
resample_bssfp(subjectID)
