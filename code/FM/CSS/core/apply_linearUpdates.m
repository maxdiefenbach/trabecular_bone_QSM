function paramsMatrix = apply_linearUpdates(paramsMatrix, rhoComplex, C_rho, C_phi)
    
    nSpecies = size(paramsMatrix, 1);
    rhoUpdate = zeros(nSpecies, 1);
    phiUpdate = zeros(nSpecies, 1);
    nRho = sum(logical(diag(C_rho)));
    nPhi = sum(logical(diag(C_phi)));
    
    for iRho = 1:nRho
        rhoUpdate = rhoUpdate + C_rho(:, iRho) .* abs(rhoComplex(iRho));
    end 
    for iPhi = 1:nPhi
        phiUpdate = phiUpdate + C_phi(:, iPhi) .* angle(rhoComplex(iPhi));
    end
    
    paramsMatrix(:, 1) = rhoUpdate;
    paramsMatrix(:, 2) = phiUpdate;

end