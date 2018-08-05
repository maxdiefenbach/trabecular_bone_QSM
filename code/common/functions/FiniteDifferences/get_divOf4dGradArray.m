function div = get_divOf4dGradArray(grad4d, voxelSize_mm, finiteDifferencesFunction)
    
    if nargin < 3
        get_finiteDifferences = @get_backwardFiniteDifferences;
    else
        get_finiteDifferences = finiteDifferencesFunction;
    end
    
    Gi = grad4d(:, :, :, 1);
    Gj = grad4d(:, :, :, 2);
    Gk = grad4d(:, :, :, 3);

    GGi = get_finiteDifferences(Gi, 1) / voxelSize_mm(1);
    GGj = get_finiteDifferences(Gj, 2) / voxelSize_mm(2);
    GGk = get_finiteDifferences(Gk, 3) / voxelSize_mm(3);
    
    div = GGi + GGj + GGk;
        
end


