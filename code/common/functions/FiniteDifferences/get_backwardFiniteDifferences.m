function G = get_backwardFiniteDifferences(array, dimension)
% G = get_forwardFiniteDifferences(array, dimension, order)
% with the Dirichlet Boundary Condition

    assert(ndims(array) == 3);
    
    switch dimension
      case 1
        arrayShift = cat(dimension, array(1, :, :), array(1:(end-1), :, :));
        
      case 2
        arrayShift = cat(dimension, array(:, 1, :), array(:, 1:(end-1), :));

      case 3
        arrayShift = cat(dimension, array(:, :, 1), array(:, :, 1:(end-1)));
    
    end
    
    G = arrayShift - array;
    
end


