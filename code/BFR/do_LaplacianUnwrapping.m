function phaseUnwrapped_rad = do_LaplacianUnwrapping(phaseWrapped_rad, voxelSize_mm, doZeropadding)
% phaseUnwrapped_rad = do_LaplacianUnwrapping(phaseWrapped_rad, voxelSize_mm, doZeropadding)
% 
% Schofield and Zhu, 2003 M.A. Schofield and Y. Zhu, Fast phase unwrapping
% algorithm for interferometric applications, Opt. Lett. 28 (2003), pp. 1194?196
    
    if nargin < 3
        doZeropadding = false
    end
    if nargin < 2
        voxelSize_mm = [1, 1, 1]
    end
    
    nDims = ndims(phaseWrapped_rad);
    if nDims == 3
        phaseUnwrapped_rad = do_LaplacianUnwrapping3D(phaseWrapped_rad, voxelSize_mm, doZeropadding);
    
    elseif nDims == 4
        phaseUnwrapped_rad = zeros(size(phaseWrapped_rad));
        nTE = size(phaseWrapped_rad, 4);
        for iTE = 1:nTE
            phaseUnwrapped_rad(:, :, :, iTE) = do_LaplacianUnwrapping3D(phaseWrapped_rad(:, :, :, iTE), ...
                                                              voxelSize_mm, ...
                                                              doZeropadding);
        end

    end

end