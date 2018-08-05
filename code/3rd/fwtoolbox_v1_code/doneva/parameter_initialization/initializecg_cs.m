function params = initializecg_cs(sx, sy)

% Initialize parameters for ncg cs reconstruction

params.alpha =  0.01;
params.beta  =  0.6;

params.image = zeros(sx,sy);

params.TV = FD1OP();  
params.WT = Wavelet_TT('Daubechies',4,4,4);

params.TVWeight = 0.04;	    % TV penalty
params.WTWeight = 0.0;      % wavelet transform l1 penalty

params.n_outer = 8;         % default number of outer iterations
params.n = 10;

params.l1Smooth = 1e-15;    % smoothing parameter of L1 norm
params.pNorm = 1;           % type of norm to use (i.e. L1 L2 etc)

params.debug = 1;

