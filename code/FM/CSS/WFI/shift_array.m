function arrayShifted = shift_array(array, iDim, dirStr, boundaryConditionStr)
    
    if nargin < 4
        boundaryConditionStr = 'replicate'; % Dirichlet, periodic, replicate
    end
    if nargin < 3
        dirStr = 'forward';
    end

    switch dirStr
      case 'forward'
        arrayShifted = get_subArray(array, iDim, 2:size(array, iDim));
        switch boundaryConditionStr
          case 'periodic'
            array2fill = get_subArray(array, iDim, 1);
          case 'replicate'
            array2fill = get_subArray(array, iDim, size(array, iDim));
          case 'Dirichlet'
            array2fill = zeros(size(get_subArray(array, iDim, 1)));
        end
        arrayShifted = cat(iDim, arrayShifted, array2fill);
      
      case 'backward'
        arrayShifted = get_subArray(array, iDim, 1:(size(array, iDim)-1));
        array2fill = get_subArray(array, iDim, size(array, iDim));
        switch boundaryConditionStr
          case 'periodic'
            array2fill = get_subArray(array, iDim, size(array, iDim));
          case 'replicate'
            array2fill = get_subArray(array, iDim, 1);
          case 'Dirichlet'
            array2fill = zeros(size(get_subArray(array, iDim, 1)));
        end
        arrayShifted = cat(iDim, array2fill, arrayShifted);
    end

end