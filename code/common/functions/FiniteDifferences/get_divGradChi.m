function divGradChi = get_divGradChi(chimap_ppm, voxelSize_mm, weighting, finiteDifferencesFunction)
    
    if nargin < 4
        get_finiteDifferences = @get_backwardFiniteDifferences;
    else
        get_finiteDifferences = finiteDifferencesFunction;
    end
    
    if nargin < 3
        weighting = 1.;
    end

    G = get_gradChi(chimap_ppm, voxelSize_mm);
    
    WGi = G(:, :, :, 1);
    WGj = G(:, :, :, 2);
    WGk = G(:, :, :, 3);

    GWGi = get_finiteDifferences(WGi, 1) / voxelSize_mm(1);
    GWGj = get_finiteDifferences(WGj, 2) / voxelSize_mm(2);
    GWGk = get_finiteDifferences(WGk, 3) / voxelSize_mm(3);
    
    divGradChi = weighting .* (GWGi + GWGj + GWGk);
        
end


