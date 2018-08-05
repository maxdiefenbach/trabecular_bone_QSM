function OutParams = get_T2s_ms(OutParams)
    
    if isfield(OutParams, 'R2s_Hz')
        T2s_ms = 1e3 ./ OutParams.R2s_Hz;
        T2s_ms(isinf(T2s_ms)) = 0;
        OutParams.T2s_ms = T2s_ms;
    end
    if isfield(OutParams, 'waterR2s_Hz')
        wT2s_ms = 1e3 ./ OutParams.waterR2s_Hz;
        wT2s_ms(isinf(wT2s_ms)) = 0;
        OutParams.waterT2s_ms = wT2s_ms;
    end
    if isfield(OutParams, 'fat_Hz')
        fT2s_ms = 1e3 ./ OutParams.fatR2s_Hz;
        fT2s_ms(isinf(fT2s_ms)) = 0;
        OutParams.fatT2s_ms = fT2s_ms;
    end

end
