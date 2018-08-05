import os
import re
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib2tikz


plt.style.use('ggplot')

def regplot(x=None, y=None, data=None, ax=None, markersize=180, txtloc=(1,1)):
    """
    regression plot of column y agains column x

    scatter with colored ROI labels
    + regresion line + CI
    + regression parameters
    + subject numbers    
    """
    df = data.copy()
    if ax is None:
        fig, ax = plt.subplots()
    colordict={1:'red', 2:'blue', 3:'green'}
    colors = df['label'].apply(lambda i: colordict[i])
    ax.scatter(df[x], df[y], c=colors, s=markersize)  # points
    sns.regplot(x=x, y=y, data=df,
                scatter=False,
                line_kws={'color': 'k'},
                ax=ax)
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(df[x], df[y])
    ax.text(*txtloc,            # regression parameters
            ('y={:.4f}x+{:.4f}\n'
             'pearson r={:.3f}\n'
             'R-squared={:.3f}\n'
             'pvalue={:.4f}').format(slope, intercept, rvalue, rvalue**2, pvalue),
            horizontalalignment='right',
            verticalalignment='top',
            Transform=ax.transAxes,
            fontsize=14)
    for i in df.index:          # subject numbers
        ax.text(df[x][i], df[y][i], df['subject'][i],
                color='white',
                horizontalalignment='center', verticalalignment='center')
    return rvalue, pvalue



df = pd.read_csv('../../../../scripts/labelstats/labelstats_total.csv')

# sns.pairplot(df.query('label in [1, 3] and '
#                       'WFI == "CSS" and '
#                       'LaplacianUnwrapping == True and '
#                       'BFR == "LBV" and '
#                       'concomCorr == True'),
#              vars=['apparentBVTV', 'mean_R2s', 'mean_chi_MEDIl1TVnesta'],
#              hue='label',
#              kind='scatter',
#              diag_kind='hist')
# plt.savefig('fig_reg.pdf', bbox_inches='tight')

df = df.query('label in [1, 3] and '
              'WFI == "CSS" and '
              'LaplacianUnwrapping == True and '
              'BFR == "LBV" and '
              'concomCorr == True')
df.to_csv('figure.dat', index=False)

plt.close('all')
fig, axs = plt.subplots(2, 2, figsize=(15, 15))
axs[0][1].axis('off')
regplot(x='apparentBVTV', y='mean_R2s', data=df,
        ax=axs[0][0], txtloc=(1,0.3))
regplot(x='apparentBVTV', y='mean_chi_MEDIl1TVnesta', data=df,
        ax=axs[1][0])
regplot(x='mean_R2s', y='mean_chi_MEDIl1TVnesta', data=df,
        ax=axs[1][1])
axs[1][1].set_ylim([-1.05, 0.18])
axs[1][1].set_xlim([45, 290])
matplotlib2tikz.save('draft_regression.tex')


# plt.close('all')
fig, axs = plt.subplots(2, 2, figsize=(15, 15))
axs[0][1].axis('off')
regplot(x='apparentBVTV', y='mean_R2s', data=df,
        ax=axs[0][0], txtloc=(1,0.3))
regplot(x='apparentBVTV', y='mean_chi_closedFormL2', data=df,
        ax=axs[1][0])
regplot(x='mean_R2s', y='mean_chi_closedFormL2', data=df,
        ax=axs[1][1])
matplotlib2tikz.save('draft_regression_closedFormL2.tex')


# plt.close('all')
fig, axs = plt.subplots(2, 2, figsize=(15, 15))
axs[0][1].axis('off')
regplot(x='apparentBVTV', y='mean_R2s', data=df,
        ax=axs[0][0], txtloc=(1,0.3))
regplot(x='apparentBVTV', y='mean_chi_MEDIl2TVcg', data=df,
        ax=axs[1][0])
regplot(x='mean_R2s', y='mean_chi_MEDIl2TVcg', data=df,
        ax=axs[1][1])
matplotlib2tikz.save('draft_regression_MEDIl2TVcg.tex')


# plt.close('all')
fig, axs = plt.subplots(2, 2, figsize=(15, 15))
dfno10 = df.query('subject != 10')
axs[0][1].axis('off')
regplot(x='apparentBVTV', y='mean_R2s', data=dfno10,
        ax=axs[0][0], txtloc=(1,0.3))
regplot(x='apparentBVTV', y='mean_chi_MEDIl2TVcg', data=dfno10,
        ax=axs[1][0])
regplot(x='mean_R2s', y='mean_chi_MEDIl2TVcg', data=dfno10,
        ax=axs[1][1])
matplotlib2tikz.save('draft_regression_MEDIl2TVcg_outlierfree.tex')


# import SimpleITK as sitk
# import os
# os.environ["SITK_SHOW_COMMAND"] = "/Applications/ITK-SNAP.app/Contents/MacOS/ITK-SNAP"

# # Image = sitk.ReadImage('../../data/reconstructed/subject01/20161104_161845_1302_20161104calcaneusexporttool4.2_WIPbFFEhiresSENSE_14.nii.gz')
# Image = sitk.ReadImage('../../data/reconstructed/subject01/20161104_161845_1302_echoMIP.nii.gz')
# sitk.Show(Image)

# Label = sitk.ReadImage(flabel)
    
# Overlay = sitk.LabelOverlay(sitk.Cast(sitk.RescaleIntensity(Image), Label.GetPixelID()),
#                             Label,
#                             opacity=0.3)

