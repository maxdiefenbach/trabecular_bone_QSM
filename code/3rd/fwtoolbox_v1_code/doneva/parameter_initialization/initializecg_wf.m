function     param = initializecg_wf(imsize, f_wf, rel_amp, t, mask)


% parameters for the backtracking line search

param.alpha         = 0.01;
param.beta          = 0.5;
param.maxIterLS     = 100;


param.pnorm_im      = 1;
param.pnorm_off     = 2;
param.maxIter       = 10;
param.threshold     = 0;

param.imsize        = imsize;
param.f_wf          = f_wf;     
param.rel_amp       = rel_amp;
param.t             = t;


param.FD2           = FD2OP();
param.FD1           = FD1OP();
param.WT            = Wavelet_TT('Daubechies',4,4,4);
param.dataWeight    = 1; 



param.FD1Weight     = 0.04;
param.FD2Weight     = 0.0000002;
param.WTWeight      = 0.0;


param.mu            = 0.000001;
param.TVweight      = 0.01;

param.debug         = 0;
param.display       = 0;
param.mask          = mask;

param.alpha_n       = 0;

% parameters for the nonlinear iterations
param.num_inner_iter_wfcs = 5;
param.num_outer_iter_wfcs = 10;

