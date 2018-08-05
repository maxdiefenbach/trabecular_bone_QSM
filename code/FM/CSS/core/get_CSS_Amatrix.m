function A = get_CSS_Amatrix(TE_s, paramsMatrix)
% A = get_CSS_Amatrix(TE_s, paramsMatrix)
% 
% without phases P = diag(exp(1j * phi_1), ..., exp(1j * phi_M))

    nTE = length(TE_s);
    nSpecies = size(paramsMatrix, 1);
    A = complex(zeros(nTE, nSpecies));
    for iSpecies = 1:nSpecies
        freq_Hz = paramsMatrix(iSpecies, 3);
        R2s_Hz = paramsMatrix(iSpecies, 4);
        for iTE = 1:nTE
            t = TE_s(iTE); 
            A(iTE, iSpecies) = exp(2j * pi * freq_Hz * t - R2s_Hz * t);
        end
    end

end