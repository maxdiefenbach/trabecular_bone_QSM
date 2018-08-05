import numpy as np
import SimpleITK as sitk
import pandas as pd
from skimage.util.montage import montage2d
import matplotlib.pyplot as plt
import scipy.stats as stats
from pprint import pprint


def montage_ndimage(ndimg, **kwargs):
    """create montage mpl plot for 3D ndimage
    in gray (ndim == 3)
    or color (ndim == 4)
    """
    if np.ndim(ndimg) == 3:
        m = montage2d(ndimg, fill=0)
    elif np.ndim(ndimg) == 4:  # + rgb dimension
        m = np.stack((montage2d(np.squeeze(ndimg[..., 0]), fill=0),
                      montage2d(np.squeeze(ndimg[..., 1]), fill=0),
                      montage2d(np.squeeze(ndimg[..., 2]), fill=0)),
                     axis=2)
    im = plt.imshow(m, **kwargs)
    return im


def get_BMDcalibration(fct, flabel):
    '''
    implementation of equation from manual of:
    QCT PRO™ Bone Minearal Densitometry Software Phantom Module
    Version 4.0 - Revision 20050430
    Copyright 2005 Mindways Software, Inc. All rights reserved.
    
    This shitty manual has very inconsitent equations due to mismatch in physical units.
    | Reference Rod | H2O density (mg/cc) | H2O stdev | K2HPO4 density (mg/cc) | K2HPO4 stdev |
    | A             |              1012.2 |       2.3 |                  -51.8 |          0.1 |
    | B             |              1057.0 |       1.9 |                  -53.4 |          0.1 |
    | C             |              1103.6 |       1.7 |                   58.9 |          0.1 |
    | D             |              1119.5 |       1.8 |                  157.0 |          0.3 |
    | E             |               923.2 |       2.1 |                  375.8 |          0.9 |
    
    (A) (B) (C) (D) (E)
    image orientation in radiological convention from inferior -> suprior
    left -> right
    dark -> bright
    (A) is the darkest, lowest HU
    (E) is the brightest, highest HU
    
    mu_ROI_phantom = rho_H2O + sigma_ref rho_K2HPO4 + beta_ref
    sigma_CT = sigma_ref - 0.2174
    beta_CT = beta_ref + 999.6
    rho_K2HPO4 = (mu_ROI - beta_CT) / sigma_CT
    '''
    df = get_labelstats_df(fct, flabel)
    df = df[df.label != 0].sort_values(by='mean')
    df['Reference Rod'] = ['A', 'B', 'C', 'D', 'E']
    df = df.set_index('Reference Rod')

    df_phantom = pd.DataFrame({
        'Reference Rod': ['A', 'B', 'C', 'D', 'E'],
        'H2O density (mg/cc)': [1012.2, 1057.0, 1103.6, 1119.5, 923.2],
        'H2O stdev': [2.3, 1.9, 1.7, 1.8, 2.1],
        'K2HPO4 density (mg/cc)': [-5.18, -53.4, 58.9, 157.0, 375.8],
        'K2HPO4 stdev': [0.1, 0.1, 0.1, 0.3, 0.9]})
    df_phantom = df_phantom.set_index('Reference Rod')

    y = df['mean'] - df_phantom['H2O density (mg/cc)']
    x = df_phantom['K2HPO4 density (mg/cc)']
    fit = stats.linregress(x, y)

    sigma_ref = fit.slope
    beta_ref = fit.intercept

    sigma_CT = sigma_ref - 0.2174
    beta_CT = beta_ref + 999.6

    print('sigma_CT = {}, beta_CT = {}'.format(sigma_CT, beta_CT))
    return sigma_CT, beta_CT


def get_labelstats_df_list(fimage_list, flabel_list):
    """loop over lists of image and label files and
    extract label statisics as pandas.DataFrame
    """
    if np.ndim(fimage_list) == 0:
        fimage_list = [fimage_list]
    if np.ndim(flabel_list) == 0:
        flabel_list = [flabel_list]

    columns = ['imagefile', 'labelfile', 'label', 'mean', 'var', 'min', 'max',
               'median', 'count', 'sum', 'boundingbox', 'voxels']
    DF = pd.DataFrame(columns=columns)
    for fimage in fimage_list:
        for flabel in flabel_list:
            df = get_labelstats_df(fimage, flabel)
            df['imagefile'] = fimage
            df['labelfile'] = flabel
        DF = DF.append(df)
    return DF


def get_labelstats_df(image, label):
    """extract image label statisics
    input:  image -- sitkImage or filename
            label -- sitkImage or filename
    output: df    -- pandas.DataFrame
    """
    labelstatsFilter = sitk.LabelStatisticsImageFilter()
    labelstatsFilter.SetUseHistograms(True)

    image = get_sitkImage(image)
    label = get_sitkImage(label)

    if label.GetDirection() != image.GetDirection() or \
       label.GetOrigin() != image.GetOrigin() or \
       label.GetSpacing() != image.GetSpacing():
        label.CopyInformation(image)
        print('Warning: copied information from image to label.')

    labelstatsFilter.Execute(image, label)

    df = pd.DataFrame(columns=['label', 'mean', 'var', 'min', 'max',
                               'median', 'count', 'sum', 'boundingbox',
                               'voxels'])
    for l in labelstatsFilter.GetLabels():
        df.loc[l] = [l,
                     labelstatsFilter.GetMean(l),
                     labelstatsFilter.GetVariance(l),
                     labelstatsFilter.GetMinimum(l),
                     labelstatsFilter.GetMaximum(l),
                     labelstatsFilter.GetMedian(l),
                     labelstatsFilter.GetCount(l),
                     labelstatsFilter.GetSum(l),
                     labelstatsFilter.GetBoundingBox(l),
                     sitk.GetArrayFromImage(image)\
                     [sitk.GetArrayFromImage(label) == l]]
    return df


def get_sitkImage(img_or_str):
    if isinstance(img_or_str, sitk.Image):
        return img_or_str
    elif isinstance(img_or_str, str):
        return sitk.ReadImage(img_or_str)
