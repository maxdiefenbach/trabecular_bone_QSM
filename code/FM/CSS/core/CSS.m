function [paramsMatrixArray, residualArray, iterArray, precisionArray, nSkippedVoxels] = ...
        CSS(signalArray, TE_s, paramsMatrixArray, constrainMatrices, tol, iterMax, skip_threshold_percent) %#codegen
% [paramsMatrixArray, residualArray, iterArray, precisionArray, nSkippedVoxels] = ...
%    CSS(signalArray, TE_s, paramsMatrixArray, constrainMatrices, tol, iterMax, skip_threshold_percent)
%
% Implementation of the VARPRO solver for chemical species separation (CSS)
% 
% described in 
% Diefenbach, M. N., Ruschke, S., & Karampinos, D. C., 
% A generalized formulation for parameter estimation in mr signals of multiple chemical species, 
% In , Proceedings 25. Annual Meeting International Society for Magnetic Resonance in Medicine 
% (pp. 5181) (2017). Honolulu, Hawaii: \url{http://dev.ismrm.org/2017/5181.html}.
% 
% Input:
%         signalArray            -- nVoxel x nTE
%         TE_s                   -- 1 x nTE
%         paramsMatrixArray      -- nVoxel x nSpecies x 4 for initialization
%                                   for each voxel
%                                   column 1: concentrations, rho's [a.u.]
%                                   column 2: transverse angle, phi's [rad]
%                                   column 3: resonance freq offset from center freq, omega's [Hz]
%                                   column 4: transverse relaxation rates, r's [s^-1]
%                                   fill column 3 for fieldmap + chemical shift
%         constrainMatrices      -- nSpecies x nSpecies x 4
%                                   dim 3 for 1:4 is Crho, Cphi, Comega, Cr
%         tol                    -- tolerance for conversion criterium
%         iterMax                -- maximum number of iteration
%         skip_threshold_percent -- skip voxel if MIP_TE smaller than ?% of maximum
%         
% Output:
%         paramsMatrixArray      -- nVoxel x nSpecies x 4, results parameter fit
%         residualArray          -- nVoxel x 1, residual norms
%         iterArray              -- nVoxel x 1, number of iterations
%         precisionArray         -- nVoxel x 1, achieved precision
%         nSkippedVoxels         -- int, number of voxels with MIP_TE value below skip_threshold_percent
%         
% function meant for code generation
% >> mex -setup
% >> codegen CSS -o CSS_mex
    
    assert(isa(signalArray, 'double') && ...
           ~isreal(signalArray) && ...
           all(size(signalArray) >= [1, 1]) && ...
           all(size(signalArray) <= [1e8, 50]));
    assert(isa(TE_s, 'double') && ...
           isreal(TE_s) && ...
           all(size(TE_s) >= [1, 2]) && ...
           all(size(TE_s) <= [1, 50]));
    assert(isa(paramsMatrixArray, 'double') && ...
           isreal(paramsMatrixArray) && ...
           all(size(paramsMatrixArray) >= [1, 1, 4]) && ...
           all(size(paramsMatrixArray) <= [1e8, 12, 4]));
    assert(isa(constrainMatrices, 'double') && ...
           all(size(constrainMatrices) >= [1, 1, 4]) && ...
           all(size(constrainMatrices) <= [12, 12, 4]));
    assert(isa(tol, 'double') && ...
           isreal(tol) && ...
           isscalar(tol));
    assert(isa(iterMax, 'double') && ...
           isreal(iterMax) && ...
           isscalar(iterMax));
    assert(isa(skip_threshold_percent, 'double') && ...
           isreal(skip_threshold_percent) && ...
           isscalar(skip_threshold_percent));
    
    nVoxels = size(signalArray, 1);
    residualArray = zeros(nVoxels, 1);
    iterArray = zeros(nVoxels, 1);
    precisionArray = zeros(nVoxels, 1);
    residual = 0.0;
    iter = 0;
    precision = 0.0;
    signalMax = max(sqrt(sum(abs(signalArray(:)).^2, 2)));
    threshold = skip_threshold_percent/100 * signalMax;
    nSkippedVoxels = 0;
    for iVoxel = 1:nVoxels
        signal = squeeze(signalArray(iVoxel, :));
        signal = signal(:);
        if mod(iVoxel, round(nVoxels/10)) == 0
            fprintf('=');
        end
        if norm(signal) < threshold
            nSkippedVoxels = nSkippedVoxels + 1;
            continue
        end
        paramsMatrix = squeeze(paramsMatrixArray(iVoxel, :, :));
        
        [paramsMatrix, residual, iter, precision] = ...
            CSSvoxel(signal, TE_s, paramsMatrix, constrainMatrices, tol, iterMax);
        
        paramsMatrixArray(iVoxel, :, :) = paramsMatrix;
        residualArray(iVoxel) = residual;
        iterArray(iVoxel) = iter;
        precisionArray(iVoxel) = precision;
    end
    
    fprintf('\n');

end