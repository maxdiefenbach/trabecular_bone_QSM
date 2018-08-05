function jacobian = compute_multiSpeciesJacobian(TE_s, paramsMatrix, constrainMatrices)

    nTE = length(TE_s);
    nSpecies = size(paramsMatrix, 1);
    
    rho = paramsMatrix(:, 1);
    P = exp(1j * paramsMatrix(:, 2));
    rhoComplex =  rho .* P;
    
    A = get_CSS_Amatrix(TE_s, paramsMatrix);
    T = diag(TE_s);
    TA = T * A;
    C_rho = squeeze(constrainMatrices(:, :, 1));
    C_phi = squeeze(constrainMatrices(:, :, 2));
    C_omega = squeeze(constrainMatrices(:, :, 3));
    C_r = squeeze(constrainMatrices(:, :, 4));
    
    % compute derivatives
    j_rho = A * ((P * diag(C_rho).') .* C_rho);
    j_phi = 1j * A * ((rhoComplex * diag(C_phi).') .* C_phi);
    j_omega = 2j * pi * TA * ((rhoComplex * diag(C_omega).') .* C_omega);
    j_r = -TA * ((rhoComplex * diag(C_r).') .* C_r);
    
    jacobian = cat(2, j_rho, j_phi, j_omega, j_r);
end