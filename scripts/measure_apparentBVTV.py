import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt
from skimage.util.montage import montage2d
from pprint import pprint


def measure_BVTV(flabel, DF):
    fileID = os.path.basename(flabel)[0:20]
    fbssfp = os.path.join(os.path.dirname(flabel), fileID + '_bssfp_N4biasCorr.nii.gz')
    subnum = os.path.basename(os.path.dirname(flabel))

    print(subnum, fileID)

    Image_bssfp = sitk.ReadImage(fbssfp)
    Label = sitk.ReadImage(flabel)

    mrimg = sitk.GetArrayFromImage(Image_bssfp)
    label = sitk.GetArrayFromImage(Label)

    fbonemask = fbssfp.replace('.nii.gz', '_Segmentation-label.nrrd')
    LabelBone = sitk.ReadImage(fbonemask)
    maskBone = sitk.GetArrayFromImage(LabelBone) > 0

    label_BVTV = []
    for l in range(1, 4):
        maskROI = label == l
        BVTV = get_apparentBVTV(mrimg, maskROI, maskBone)
        label_BVTV.append([l, BVTV])

    df = pd.DataFrame(data=label_BVTV, columns=['label', 'apparentBVTV'])
    df['subjectID'] = subnum
    df['fileID'] = fileID
    df = df[columns]
    print(df)
    df.to_csv(fbssfp.replace('bssfp_N4biasCorr.nii.gz', 'BVTV.csv'), index=False)
    DF = DF.append(df)
    return DF
    


def get_apparentBVTV(mrimg, maskROI, maskBone, **kwargs):
    ''' compute bone volume to total volume ratio
    
    Majumdar, S., Genant, H. K., Grampp, S., Newitt, D. C., 
    Truong, V., Lin, J. C., & Mathur, A., 
    Correlation of trabecular bone structure with age, bone mineral density, 
    and osteoporotic status: in vivo studies in the distal radius using high
    resolution magnetic resonance imaging, 
    Journal of Bone and Mineral Research, 12(1), 111–118 (1997).  
    http://dx.doi.org/10.1359/jbmr.1997.12.1.111
    '''
    inverted = (np.max(mrimg) - mrimg)
    nbins = kwargs.get('nbins', 'auto')
    hist, binEdges = np.histogram(inverted[maskROI], bins=nbins)
    inds2crop = kwargs.get('inds2crop', 0) + 1
    indMaxCount = np.argmax(hist[:-inds2crop])
    indHalfMaxCount = np.argmin(np.abs(hist[:(indMaxCount+1)] - hist[indMaxCount] / 2))
    I_L = binEdges[indHalfMaxCount+1]
    I_b = np.mean(inverted[maskBone])
    I_r = np.mean(inverted[maskROI])
    BVTV = (I_r - I_L) / (I_b - I_L)
    if kwargs.get('verbose', False):
        print("nbins = {}\ninds2crop = {}".format(nbins, inds2crop))
        print("Calculate bone volume to total volume ratio\nBV/TV = {}\nby BV/TV = (I_r - I_L) / (I_b - I_L) with\nI_L = {}\nI_b = {}\nI_r = {}\n".format(BVTV, I_L, I_b, I_r))
    # assert BVTV >= 0
    if kwargs.get('showHist', False):
        plt.figure()
        plt.bar(binEdges[1:], hist)
        plt.plot(binEdges[indMaxCount+1], hist[indMaxCount], 'or')
        plt.plot(binEdges[indHalfMaxCount+1], hist[indHalfMaxCount], 'ob')
        if inds2crop != 1:
            plt.axvline(x=binEdges[len(binEdges)-inds2crop], color='black')
        plt.title('ROI histogram of inverted image')
        plt.show()
    if kwargs.get('returnHist', False):
        return(BVTV, hist, binEdges)
    else:
        return(BVTV)


datadir = '../data/reconstructed'

label_files = [os.path.join(r, f) for r, d, fs in os.walk(datadir)
               for f in fs if f.endswith('Segmentation_resamp-label.nrrd')]
pprint(label_files)

columns = ['subjectID', 'fileID', 'label', 'apparentBVTV']
DF = pd.DataFrame(columns=columns)

for flabel in label_files:
    DF = measure_BVTV(flabel, DF)

DF.to_csv('../data/analysis/BVTV.csv', index=False)
