function array_flat = flatten3d(array)
% flatten the first three array dimensions
    
    if isvector(array)
        array_flat = array;
    else
        arraySize = [size(array, 1), size(array, 2), ...
                     size(array, 3), size(array, 4)];
        % adds singelton dimenstion if ndim(array) = 2 or 4

        Length = prod(arraySize(1:3));
        arraySizeNew = [Length, arraySize(4:end)];
        if length(arraySizeNew) == 1
            arraySizeNew(2) = 1;
        end

        array_flat = reshape(array, arraySizeNew);
    end
    
end % EOF