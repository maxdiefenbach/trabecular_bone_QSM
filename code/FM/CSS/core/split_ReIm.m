function array = split_ReIm(array)
% array = split_ReIm(array)
    array = cat(1, real(array), imag(array));
end