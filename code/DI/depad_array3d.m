function array3dSmall = depad_array3d(array3dBig, padsize)
    
    assert(ndims(array3dBig) == 3)
    assert(length(padsize) == 3)

    matrixSize = size(array3dBig);

    x = (padsize(1)+1):(matrixSize(1)-padsize(1));
    y = (padsize(2)+1):(matrixSize(2)-padsize(2));
    z = (padsize(3)+1):(matrixSize(3)-padsize(3));
    
    array3dSmall = array3dBig(x, y, z);
end
