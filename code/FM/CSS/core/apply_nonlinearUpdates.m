function paramsMatrix = apply_nonlinearUpdates(paramsMatrix, updates, C_rho, C_phi, C_omega, C_r)
    
    nSpecies = size(paramsMatrix, 1);
    
    nRho = 0;
    nPhi = 0;
    % nRho = sum(logical(diag(C_rho)));
    % rhoUpdate = zeros(nSpecies, 1);
    % for iRho = 1:nRho
    %     rhoUpdate = rhoUpdate + C_rho(:, iRho) .* updates(iRho);
    % end 
    
    % nPhi = sum(logical(diag(C_phi)));
    % phiUpdate = zeros(nSpecies, 1);
    % for iPhi = 1:nPhi
    %     phiUpdate = phiUpdate + C_phi(:, iPhi) .* updates(nRho + iPhi);
    % end

    nOmega = sum(logical(diag(C_omega)));
    omegaUpdate = zeros(nSpecies, 1);
    for iOmega = 1:nOmega
        omegaUpdate = omegaUpdate + C_omega(:, iOmega) * updates(nRho + nPhi + iOmega);
    end
    
    nR = sum(logical(diag(C_r)));
    rUpdate = zeros(nSpecies, 1);
    for iR = 1:nR
        rUpdate = rUpdate + C_r(:, iR) .* updates(nRho + nPhi + nOmega + iR); 
    end

    % paramsMatrix(:, 1) = paramsMatrix(:, 1) + rhoUpdate;
    % paramsMatrix(:, 2) = paramsMatrix(:, 2) + phiUpdate;
    paramsMatrix(:, 3) = paramsMatrix(:, 3) + omegaUpdate;
    paramsMatrix(:, 4) = paramsMatrix(:, 4) + rUpdate;

end