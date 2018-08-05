function hessian = compute_multiSpeciesHessian(TE_s, paramsMatrix, constrainMatrices)

    nTE = length(TE_s);
    nSpecies = size(paramsMatrix, 1);
    
    rho = paramsMatrix(:, 1);
    drho = diag(rho);
    P = exp(1j * paramsMatrix(:, 2));
    A = get_CSS_Amatrix(TE_s, paramsMatrix);
    T = diag(TE_s);
    TA = T * A;
    AP = A * P;
    TAP = T * A * P;
    C_rho = squeeze(constrainMatrices(:, :, 1));
    C_phi = squeeze(constrainMatrices(:, :, 2));
    C_omega = squeeze(constrainMatrices(:, :, 3));
    C_r = squeeze(constrainMatrices(:, :, 4));
    
    % compute second derivatives
    h_rho_rho = zeros(size(C_rho));
    h_rho_phi = 1j * AP * (C_rho .* C_phi);
    h_rho_omega = 1j * TAP * (C_rho .* C_omega);
    h_rho_r = 1j * TAP * (C_rho .* C_r);
    h_phi_phi = -AP * drho * C_phi;
    h_phi_omega = -TAP * drho * (C_phi .* C_omega):
    h_phi_r = 1j * TAP * drho * (C_phi .* C_r);
    h_omega_omega = -T^2 * AP * drho * C_omega;
    h_omega_r = -1j * T^2 * AP * drho * (C_omega .* C_r);
    h_r_r = -T^2 * AP * drho * C_r;
    
    hessian = cat(2, h_rho_rho, h_rho_phi, h_rho_omega, h_rho_r);
end