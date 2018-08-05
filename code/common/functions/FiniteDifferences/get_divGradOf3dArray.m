function div = get_divGradOf3dArray(array3d, voxelSize_mm, finiteDifferencesFunction)

    if nargin < 3
        get_finiteDifferences = @get_forwardFiniteDifferences;
    else 
        get_finiteDifferences = finiteDifferencesFunction;
    end
    
    grad = get_gradOf3dArray(array3d, voxelSize_mm);
    
    div = get_divOf4dGradArray(grad, voxelSize_mm);
 
end
