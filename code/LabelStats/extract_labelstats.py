import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt
from pprint import pprint
from mySimpleITK import *


datadir = '../../data/'

label_files = [os.path.join(r, f) for r, d, fs in os.walk(datadir)
               for f in fs if f.endswith('echoMIP_Segmentation-label.nrrd')
               if not r.endswith('/CT')]
pprint(label_files)

columns = ['imagefile', 'labelfile', 'label', 'mean', 'var', 'min', 'max',
           'median', 'count', 'voxels']
DF = pd.DataFrame(columns=columns)

for flabel in label_files:
    fileID = os.path.basename(flabel)[0:20]
    subjectDir = os.path.dirname(flabel)
    subjectID = os.path.basename(subjectDir)

    print(flabel)
    imagefile_list = [os.path.join(r, f) for r, d, fs in os.walk(subjectDir)
                      for f in fs
                      if f.endswith('.nii.gz')
                      # if 'echoMIP' not in f
                      if 'bFFE' not in f
                      and 'bssfp' not in f
                      if not r.endswith('/CT')]
    pprint(imagefile_list)

    df = get_labelstats_df_list(imagefile_list, flabel)[columns]
    DF = DF.append(df, ignore_index=True)

    fcsv = os.path.join(subjectDir, subjectID + '_' + fileID +
                        '_labelstats.csv')
    df[columns].to_csv(fcsv, index=False)
    print('Wrote {}'.format(fcsv))

DF = DF[columns]
DF.head()
# fcsv = '../../data/analysis/labelstats.csv'
# DF.to_csv(fcsv, index=False)
fcsv = './labelstats.csv'
fhdf = fcsv.replace('csv', 'h5')
# DF.to_csv(fcsv, index=False)
# DF.to_hdf(fhdf, '/labelstats', index=False)
# print('Wrote {} \n {}'.format(fcsv, fhdf))
