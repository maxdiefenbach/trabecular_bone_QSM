function fd = finite_difference(array, iDim, dirStr, boundaryConditionStr)
    
    if nargin < 4
        boundaryConditionStr = 'replicate';
    end
    if nargin < 3
        dirStr = 'forward';
    end
    
    if isequal(dirStr, 'central')
        farr = shift_array(array, iDim, 'forward', boundaryConditionStr);
        barr = shift_array(array, iDim, 'backward', boundaryConditionStr);
        fd = 0.5 * (farr - barr);
    else
        arrayShifted = shift_array(array, iDim, dirStr, boundaryConditionStr);
        fd = arrayShifted - array;
    end

end