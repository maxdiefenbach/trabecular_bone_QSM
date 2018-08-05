import SimpleITK as sitk
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append('/Users/maxdiefenbach/programs/BMRR/Postprocessing/LabelStatistics')
from mySimpleITK import get_labelstats_df
import scipy.stats as stats


#################### patient 1 ####################
label = sitk.ReadImage('../../../../data/reconstructed/subjectCT16/CT/CT_ROIs-label.nrrd')
label.SetOrigin((0, 0, 0))
ct = sitk.ReadImage('../../../../data/reconstructed/subjectCT16/CT/20170216_Anonymous_Male_1939_anonymous_404_rightFoot_resamp2echoMIP.nrrd')
ct.SetOrigin((0, 0, 0))
chi = sitk.ReadImage('../../../../data/reconstructed/subjectCT16/20170216_113523_0302_chi_closedFormL2_LBV_unwrap_CSS.nii.gz')
chi.SetOrigin((0, 0, 0))
r2s = sitk.ReadImage('../../../../data/reconstructed/subjectCT16/20170216_113523_0302_R2s_CSS.nii.gz')
r2s.SetOrigin((0, 0, 0))

df1_ct = get_labelstats_df(ct, label)
df1_ct = df1_ct[df1_ct['label']!=0]
df1_chi = get_labelstats_df(chi, label)
df1_chi = df1_chi[df1_chi['label']!=0]
df1_r2s = get_labelstats_df(r2s, label)
df1_r2s = df1_r2s[df1_r2s['label']!=0]
df1 = pd.DataFrame({'subject': 1,
                   'label': df1_ct['label'],
                   'mean_ct': df1_ct['mean'],
                   'var_ct': df1_ct['var'],
                   'mean_r2s': df1_r2s['mean'],
                   'var_R2s_Hz': df1_r2s['var'],
                   'mean_chi': df1_chi['mean'],
                   'var_chi_Hz': df1_chi['var']})
print(df1)


#################### patient 2 ####################
label = sitk.ReadImage('../../../../data/reconstructed/subjectCT17/CT/CT_ROIs-label.nrrd')
label.SetOrigin((0, 0, 0))
ct = sitk.ReadImage('../../../../data/reconstructed/subjectCT17/CT/20170407_Anonymous_Female_1946_anonymous_205_rightFoot_resamp2echoMIP.nrrd')
ct.SetOrigin((0, 0, 0))
chi = sitk.ReadImage('../../../../data/reconstructed/subjectCT17/20170407_115911_0302_chi_closedFormL2_LBV_unwrap_CSS.nii.gz')
chi.SetOrigin((0, 0, 0))
r2s = sitk.ReadImage('../../../../data/reconstructed/subjectCT17/20170407_115911_0302_R2s_CSS.nii.gz')
r2s.SetOrigin((0, 0, 0))

df2_ct = get_labelstats_df(ct, label)
df2_ct = df2_ct[df2_ct['label']!=0]
df2_chi = get_labelstats_df(chi, label)
df2_chi = df2_chi[df2_chi['label']!=0]
df2_r2s = get_labelstats_df(r2s, label)
df2_r2s = df2_r2s[df2_r2s['label']!=0]
df2 = pd.DataFrame({'subject': 2,
                   'label': df2_ct['label'],
                   'mean_ct': df2_ct['mean'],
                   'var_ct': df2_ct['var'],
                   'mean_r2s': df2_r2s['mean'],
                   'var_R2s_Hz': df2_r2s['var'],
                   'mean_chi': df2_chi['mean'],
                   'var_chi_Hz': df2_chi['var']})
print(df2)

df = df1.append(df2)
df['ROI'] = df['label'].apply(lambda x: 1 if x <= 10 else 2)
df.reset_index(inplace=True)
df.to_csv('CT_ROIs_labelstats.csv', index=False)
df.to_csv('figure_clfl2.dat', index=False)
print(df)


plt.close('all')
plt.style.use('ggplot')

# plt.figure()
# sns.regplot('mean_ct', 'mean_r2s', df)

# plt.figure()
# sns.regplot('mean_ct', 'mean_chi', df)

# plt.figure()
# sns.regplot('mean_r2s', 'mean_chi', df)


# sns.lmplot(x='mean_ct', y='mean_chi', hue='subject', data=df)
# sns.lmplot(x='mean_ct', y='mean_r2s', hue='subject', data=df)
# sns.lmplot(x='mean_ct', y='mean_chi', hue='subject', data=df)


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
    colors = df['ROI'].apply(lambda i: colordict[i])
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
            Transform=ax.transAxes)
    for i in df.index:          # subject numbers
        ax.text(df[x][i], df[y][i], df['subject'][i],
                color='white',
                horizontalalignment='center', verticalalignment='center')
    return rvalue, pvalue

# plt.close('all')
# regplot(x='mean_ct', y='mean_chi', data=df)
# regplot(x='mean_ct', y='mean_r2s', data=df)
# regplot(x='mean_r2s', y='mean_chi', data=df)

plt.close('all')
fig, axs = plt.subplots(2, 2, figsize=(12, 12))
axs[0][1].axis('off')
regplot(x='mean_ct', y='mean_r2s', data=df,
        ax=axs[0][0], txtloc=(1,0.3))
regplot(x='mean_ct', y='mean_chi', data=df,
        ax=axs[1][0])
regplot(x='mean_r2s', y='mean_chi', data=df,
        ax=axs[1][1])
axs[1][1].set_xlim([45, 230])
axs[1][1].set_ylim([-1.25, 0.39])
import matplotlib2tikz
matplotlib2tikz.save('draft_regression_closedFormL2.tex')
