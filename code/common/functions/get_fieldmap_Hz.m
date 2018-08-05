function WFIparams = get_fieldmap_Hz(ImDataParams, Options)
    
    M = ImDataParams.signal;
    TE = ImDataParams.TE_s;
    
    [p1, dp1, relres, p0] = Fit_ppm_complex_TE(M, TE);
    %   output
    %   p1 - field map, may need further unwrapping
    %   dp1 - a priori error estimate
    %   relres - relative residual
    %   p0 - initla phase
    
    fieldmap_Hz = p1 * 1e-6 * ImDataParams.centerFreq_Hz;
    errorEstimate = dp1;
    residual = relres;
    initialPhase_rad = p0;
    
    WFIparams = pack_vars2struct(fieldmap_Hz, ...
                                 errorEstimate, ...
                                 residual, ...
                                 initialPhase_rad);

end