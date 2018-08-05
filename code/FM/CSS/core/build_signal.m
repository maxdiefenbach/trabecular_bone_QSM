function signal = build_signal(TE_s, paramsMatrix)
% signal = build_signal(TE_s, paramsMatrix)

    nTE = length(TE_s);
    nSpecies = size(paramsMatrix, 1);
    signal = complex(zeros(nTE, 1));
    for iTE = 1:nTE
        t = TE_s(iTE);
        for iSpecies = 1:nSpecies
            magnitude = paramsMatrix(iSpecies, 1);
            phase = paramsMatrix(iSpecies, 2);
            freq_Hz = paramsMatrix(iSpecies, 3);
            R2s_Hz = paramsMatrix(iSpecies, 4);
            signal(iTE) = signal(iTE) + ...
                magnitude * exp(2j * pi * freq_Hz * t - R2s_Hz * t + 1j * phase);
        end
    end
    
end