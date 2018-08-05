clear; close all; clc;

% setup
p = path;                               % for later tear down
restoredefaultpath;
addpath(genpath('../'))



filename = '../data/20180406_144018_0302_ImDataParams.mat'
I = ImDataParamsBMRR(filename)

% I = WFI_subject(filename)
I.plot_WFIparams

I.load_WFIparams('../data/20180406_144018_0302_WFIparams_CSS_unwrap.mat')

I = BFR_subject(filename)
I.plot_BFRparams

I = QSM_subject(filename)
I.plot_QSMparams


% teardown: restore previous matlab path
path(p)