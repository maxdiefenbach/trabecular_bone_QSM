function center_ijk = get_center_ijk(ImDataParams)
    
    matrixSize = get_matrixSize(ImDataParams);
    center_ijk = ceil((matrixSize + 1) ./ 2);

end