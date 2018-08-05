function grad = get_gradOf3dArray(array3d, voxelSize_mm, finiteDifferencesFunction)

    if nargin < 3
        get_finiteDifferences = @get_forwardFiniteDifferences;
    else 
        get_finiteDifferences = finiteDifferencesFunction;
    end
    if nargin < 2
        voxelSize_mm = [1, 1, 1,];
    end
    
    Gi = get_finiteDifferences(array3d, 1) / voxelSize_mm(1);
    Gj = get_finiteDifferences(array3d, 2) / voxelSize_mm(2);
    Gk = get_finiteDifferences(array3d, 3) / voxelSize_mm(3);

    grad = cat(4, Gi, Gj, Gk);

end
