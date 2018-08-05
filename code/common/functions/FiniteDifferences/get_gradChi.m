function gradChi4d = get_gradChi(chimap_ppm, voxelSize_mm, finiteDifferencesFunction)
% gradChi4d = get_gradChi(chimap_ppm, voxelSize_mm)

    if nargin < 3
        get_finiteDifferences = @get_forwardFiniteDifferences;
    else 
        get_finiteDifferences = finiteDifferencesFunction;
    end
    if nargin < 2
        voxelSize_mm = [1, 1, 1,];
    end
    
    Gi = get_finiteDifferences(chimap_ppm, 1) / voxelSize_mm(1);
    Gj = get_finiteDifferences(chimap_ppm, 2) / voxelSize_mm(2);
    Gk = get_finiteDifferences(chimap_ppm, 3) / voxelSize_mm(3);

    gradChi4d = cat(4, Gi, Gj, Gk);

end
