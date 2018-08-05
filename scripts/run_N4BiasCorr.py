import os
import string
import re
import numpy as np
import SimpleITK as sitk
import scipy.ndimage as ndi
import nibabel as nib
from pprint import pprint


os.environ["SITK_SHOW_COMMAND"] = "/Applications/ITK-SNAP.app/Contents/MacOS/ITK-SNAP"


def run_N4biasCorr(fimg, fout):
    numberFittingLevels = 4
    numberOfIterations = 3
    shrinkFactor = 2

    corrector = sitk.N4BiasFieldCorrectionImageFilter()
    corrector.SetMaximumNumberOfIterations([int(numberOfIterations)] * numberFittingLevels)
    corrector.AddCommand(sitk.sitkStartEvent, lambda: print("Start N4 bias field correction."))
    corrector.AddCommand(sitk.sitkEndEvent, lambda: print("Done."))

    Image = sitk.ReadImage(fimg)
    Image = sitk.Shrink(Image, [int(shrinkFactor)] * Image.GetDimension())  # without this it takes an awful amount of time

    maskImage = sitk.OtsuThreshold(Image, 0, 1, 200)

    inputImage = sitk.Cast(Image, sitk.sitkFloat32)

    outputImage = corrector.Execute(inputImage, maskImage)

    sitk.WriteImage(outputImage, fout)
    print('Wrote {}.\n'.format(fout))
    

datadir = '../data/reconstructed'
bssfp_images = [os.path.join(r, f) for r, d, fs in os.walk(datadir)
                for f in fs if re.match('.*bFFEhires.*\.nii\.gz', f)]
pprint(bssfp_images)

for fbssfp in bssfp_images[-2:]:
    fileID = os.path.basename(fbssfp)[0:20]
    fout = os.path.join(os.path.dirname(fbssfp), fileID + '_bssfp_N4biasCorr.nii.gz')
    run_N4biasCorr(fbssfp, fout)


    
