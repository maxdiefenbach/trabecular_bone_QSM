function OutParams = fit_monoExponentialRelaxation(ImDataParams, AlgoParams) 
    
    if nargin < 2
        AlgoParams = struct();
    end
    if isfield(AlgoParams, 'R2sRange_Hz')
        r2smin = AlgoParams.R2sRange_Hz(1);
        r2smax = AlgoParams.R2sRange_Hz(end);
    else
        r2smin = 0;
        r2smax = 1000;
    end
    if isfield(AlgoParams, 'precision')
        abs_err = AlgoParams.precision;
    else
        abs_err = 1e-4;
    end
    
    TE_s = ImDataParams.TE_s;
    nTE = length(TE_s);
    matrixSize = get_matrixSize(ImDataParams);
    magnitude = get_magnitude(ImDataParams);
    magnitudeFlat = reshape(magnitude, [prod(matrixSize), nTE]);
    
    func = @(tn, R2s, p_) exp(-R2s .* tn);

    [s, beta, R2s, gof] = ...
        least_squares_varpro_gss_cg(TE_s, magnitudeFlat, [], ...
                                    func, r2smin, r2smax, [], abs_err, ...
                                    'show progress');

    s = reshape(s, matrixSize);
    beta = beta;
    R2s = reshape(R2s, matrixSize);
    gof = reshape(gof, matrixSize);

    OutParams.s = s;
    OutParams.beta = beta;
    OutParams.R2s_Hz = R2s;
    OutParams.gof = gof;

end