function [paramsMatrix, residual, iter, precision]  = CSSvoxel(signal, TE_s, paramsMatrix, constrainMatrices, tol, iterMax)
% [paramsMatrix, residual, iter, precision]  = CSSvoxel(signal, TE_s, paramsMatrix, constrainMatrices, tol, iterMax)
% 
% Implementation of the VARPRO solver for chemical species separation (CSS) in only one voxel
% 
% described in 
% Diefenbach, M. N., Ruschke, S., & Karampinos, D. C., 
% A generalized formulation for parameter estimation in mr signals of multiple chemical species, 
% In , Proceedings 25. Annual Meeting International Society for Magnetic Resonance in Medicine 
% (pp. 5181) (2017). Honolulu, Hawaii: \url{http://dev.ismrm.org/2017/5181.html}.
% 
% Input:
%         signal                 -- 1 x nTE
%         TE_s                   -- 1 x nTE
%         paramsMatrix           -- nSpecies x 4 for initialization
%                                   column 1: concentrations, rho's [a.u.]
%                                   column 2: transverse angle, phi's [rad]
%                                   column 3: resonance freq offset from center freq, omega's [Hz]
%                                   column 4: transverse relaxation rates, r's [s^-1]
%                                   fill column 3 for fieldmap + chemical shift
%         constrainMatrices      -- nSpecies x nSpecies x 4
%                                   dim 3 for 1:4 is Crho, Cphi, Comega, Cr
%         tol                    -- tolerance for conversion criterium
%         iterMax                -- maximum number of iteration
%         
% Output:
%         paramsMatrix           -- nSpecies x 4, results parameter fit
%         residual               -- double, residual norms
%         iter                   -- single, number of iterations
%         precision              -- double, achieved precision
    
    nTE = length(TE_s);
    nSpecies = size(paramsMatrix, 1);

    C_rho = squeeze(constrainMatrices(:, :, 1));
    C_phi = squeeze(constrainMatrices(:, :, 2));
    C_omega = squeeze(constrainMatrices(:, :, 3));
    C_r = squeeze(constrainMatrices(:, :, 4));
    unknownParamsVec = logical(cat(1, diag(C_rho), diag(C_phi), diag(C_omega), diag(C_r)));
    unknownParamsMask = logical(reshape(unknownParamsVec, [nSpecies, 4]));
    unknownParamsVec = logical(cat(1, zeros(2*nSpecies, 1), diag(C_omega), diag(C_r)));
    unknowRhoMask = logical(diag(C_rho));

    precision = 2 * tol;
    iter = 0;
    residual = 1.0e3;
    while precision > tol & iter < iterMax
        
        % update linear params
        A = get_CSS_Amatrix(TE_s, paramsMatrix);
        % P = diag(squeeze(exp(1j .* paramsMatrix(:, 2))));
        Ar = A * C_rho(:, unknowRhoMask);
        if sum(isnan(Ar(:))) | sum(isinf(Ar(:)))
            disp('A matrix has nans or infs.')
            break
        end
        rhoComplex = Ar \ signal;
        paramsMatrix = apply_linearUpdates(paramsMatrix, rhoComplex, C_rho, C_phi);
        
        % update nonlinear params
        s = build_signal(TE_s, paramsMatrix);
        ds = signal - s;
        J = compute_multiSpeciesJacobian(TE_s, paramsMatrix, constrainMatrices);
        J = J(:, unknownParamsVec);
        updates = split_ReIm(J) \ split_ReIm(ds);
        paramsMatrix = apply_nonlinearUpdates(paramsMatrix, updates, C_rho, C_phi, C_omega, C_r);
        
        residual = norm(signal - build_signal(TE_s, paramsMatrix));
        precision = norm(updates) ./ norm(paramsMatrix(unknownParamsMask));
        iter = iter + 1;
    end

end