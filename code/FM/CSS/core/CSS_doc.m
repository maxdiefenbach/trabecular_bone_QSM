function [paramsMatrixArray, residualArray, iterArray, precisionArray, nSkippedVoxels] = ...
    CSS(signalArray, TE_s, paramsMatrixArray, constrainMatrices, tol, iterMax, skip_threshold_percent)
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
