function G = get_centralFiniteDifferences(field, dimension, order)
% G = get_centralFiniteDifferences(field, dimension, order)

    assert(ndims(field) == 3);
    
    if nargin < 3
        order = 1;
    end

    finiteDifferences = diff(field, order, dimension);
    
    switch dimension
      case 1
        first = cat(dimension, finiteDifferences(1:order, :, :), finiteDifferences);
        last  = cat(dimension, finiteDifferences, finiteDifferences((end-(order-1)):end, :, :));

      case 2
        first = cat(dimension, finiteDifferences(:, 1:order, :), finiteDifferences);
        last  = cat(dimension, finiteDifferences, finiteDifferences(:, (end-(order-1)):end, :));

      case 3
        first = cat(dimension, finiteDifferences(:, :, 1:order), finiteDifferences);
        last  = cat(dimension, finiteDifferences, finiteDifferences(:, :, (end-(order-1)):end));
    end
    
    G = (first + last) / 2;

end


