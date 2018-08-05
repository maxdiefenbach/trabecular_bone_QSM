import os
import sys
import re
import SimpleITK as sitk
from pprint import pprint


os.environ["SITK_SHOW_COMMAND"] = "/Applications/ITK-SNAP.app/Contents/MacOS/ITK-SNAP"


def run_registration(ffix, fmov, ftrafo=None):
    print('\nRegistration of\n{}\nto\n{}.'.format(fmov, ffix))

    fixed_image = sitk.ReadImage(ffix, sitk.sitkFloat32)
    moving_image = sitk.ReadImage(fmov, sitk.sitkFloat32)

    initialTransform = sitk.Euler3DTransform()
    initialTransform = sitk.CenteredTransformInitializer(fixed_image,
                                                         moving_image,
                                                         initialTransform,
                                                         sitk.CenteredTransformInitializerFilter.MOMENTS)

    registration = sitk.ImageRegistrationMethod()
    registration.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration.SetMetricSamplingStrategy(registration.RANDOM)
    registration.SetMetricSamplingPercentage(0.01)
    registration.SetInterpolator(sitk.sitkLinear)
    registration.SetOptimizerAsGradientDescent(learningRate=1.0,
                                               numberOfIterations=100,
                                               convergenceMinimumValue=1e-6,
                                               convergenceWindowSize=10)
    registration.SetOptimizerScalesFromPhysicalShift()
    registration.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
    registration.SetSmoothingSigmasPerLevel(smoothingSigmas = [2,1,0])
    registration.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    registration.SetInitialTransform(initialTransform)
    registration.AddCommand(sitk.sitkStartEvent, lambda: print("Start."))
    registration.AddCommand(sitk.sitkEndEvent, lambda: print("Done"))
    transform = registration.Execute(fixed_image, moving_image)

    if ftrafo is not None:
        transform.WriteTransform(ftrafo)
        print('Wrote {}.'.format(ftrafo))
    return transform


def run_resampling(ffix, fmov, ftrafo, fresamp=None):
    print('Resample {}\nwith trafo\n{}\nto fixed image\n{}.'.format(fmov, ftrafo, ffix))
    
    fixed_image = sitk.ReadImage(ffix, sitk.sitkFloat32)
    moving_image = sitk.ReadImage(fmov, sitk.sitkFloat32)
    transform = sitk.ReadTransform(ftrafo)

    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(fixed_image)
    resample.SetInterpolator(sitk.sitkBSpline)
    resample.SetTransform(transform)
    resample.AddCommand(sitk.sitkProgressEvent, lambda: print("\rProgress: {0:03.1f}%...".format(100*resample.GetProgress()),end=''))
    resample.AddCommand(sitk.sitkProgressEvent, lambda: sys.stdout.flush())
    resample.AddCommand(sitk.sitkEndEvent, lambda: print("Done"))
    output_image = resample.Execute(moving_image)

    if fresamp is not None:
        sitk.WriteImage(output_image, fresamp)
        print('Wrote {}.'.format(fresamp))
    return output_image


def save_checkerboard(fimg1, fimg2, fout):
    Image1 = sitk.ReadImage(fimg1, sitk.sitkFloat32)
    Image2 = sitk.ReadImage(fimg2, sitk.sitkFloat32)
    CheckerBoard = sitk.CheckerBoard(Image1, Image2, checkerPattern = [8, 8, 8])
    sitk.WriteImage(CheckerBoard[:, :, cb.GetDepth()//2], fout)
    print('Wrote {}.'.format(fout))
    return CheckerBoard


datadir = '../data/reconstructed'
fixed_images = [os.path.join(r, f) for r, d, fs in os.walk(datadir)
                for f in fs if re.match('.*bFFE.*hires.*\.nii\.gz', f)]
pprint(fixed_images)

for ffix in fixed_images[-2:]:
    fileID = os.path.basename(ffix)[0:20]
    if fileID == '20170523_172947_1102':
        continue
    fmov = os.path.join(os.path.dirname(ffix), fileID + '_echoMIP.nii.gz')
    ftrafo = os.path.join(os.path.dirname(ffix), fileID + '_transform_mip2bssfp.tfm')
    run_registration(ffix, fmov, ftrafo)



