function FOV_mm = get_FOV_mm(ImDataParams)
% FOV_mm = get_FOV_mm(ImDataParams)

matrixSize = get_matrixSize(ImDataParams);
voxelSize_mm = ImDataParams.voxelSize_mm;
FOV_mm = matrixSize .* voxelSize_mm;

end