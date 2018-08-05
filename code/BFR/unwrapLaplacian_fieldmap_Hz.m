function fieldmapUnwrapped_Hz = unwrapLaplacian_fieldmap_Hz(ImDataParams, fieldmap_Hz, doZeropadding)
    
    if nargin < 3
        doZeropadding = false
    end

    matrixSize = size(fieldmap_Hz);
    voxelSize_mm = ImDataParams.voxelSize_mm;
    TE_s = ImDataParams.TE_s;
    nTE = length(TE_s);
    
    fieldmapUnwrapped_Hz = zeros([matrixSize, nTE]);
    for iTE = 1:nTE
        t_s = TE_s(iTE);
        phaseWrapped_rad = 2 * pi * fieldmap_Hz * t_s;
        phaseUnwrapped_rad = do_LaplacianUnwrapping3D(phaseWrapped_rad, ...
                                                      voxelSize_mm, ...
                                                      doZeropadding);
        fieldmapUnwrapped_Hz(:, :, :, iTE) = phaseUnwrapped_rad ./ (2 * pi * t_s);
    end

end