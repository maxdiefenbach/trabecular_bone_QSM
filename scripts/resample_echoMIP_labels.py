import os
import SimpleITK as sitk
import numpy as np
from skimage.util.montage import montage2d
import matplotlib.pyplot as plt
from pprint import pprint


def resample_label(flabel):
    fileID = os.path.basename(flabel)[0:20]
    print(fileID)
    fbssfp = os.path.join(os.path.dirname(flabel), fileID + '_bssfp_N4biasCorr.nii.gz')
    fmip = os.path.join(os.path.dirname(flabel), fileID + '_echoMIP.nii.gz')
    ftrafo = os.path.join(os.path.dirname(flabel), fileID + '_transform_mip2bssfp.tfm')

    Image_bssfp = sitk.ReadImage(fbssfp)
    Image_mip = sitk.ReadImage(fmip)
    Label = sitk.ReadImage(flabel)
    transform = sitk.ReadTransform(ftrafo)
    
    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(Image_bssfp)
    resample.SetTransform(transform)
    resample.SetInterpolator(sitk.sitkNearestNeighbor)
    Label_resamp = resample.Execute(Label)

    fout = flabel.replace('Segmentation-label', 'Segmentation_resamp-label')
    sitk.WriteImage(Label_resamp, fout)
    print('Wrote {}.'.format(fout))


datadir = '../data/reconstructed'

label_files = [os.path.join(r, f) for r, d, fs in os.walk(datadir)
               for f in fs
               if f.endswith('Segmentation-label.nrrd')]
pprint(label_files)

for flabel in label_files[-2:]:
    resample_label(flabel)
