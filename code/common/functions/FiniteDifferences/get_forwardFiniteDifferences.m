function G = get_forwardFiniteDifferences(array3d, dimension)
% G = get_forwardFiniteDifferences(array3d, dimension, order)

    assert(ndims(array3d) == 3);
    
    switch dimension
      case 1
        fieldShifted = cat(dimension, array3d(2:end, :, :), array3d(end, :, :));

      case 2
        fieldShifted = cat(dimension, array3d(:, 2:end, :), array3d(:, end, :));

      case 3
        fieldShifted = cat(dimension, array3d(:, :, 2:end), array3d(:, :, end));
    end
    
    G = fieldShifted - array3d;
    
end