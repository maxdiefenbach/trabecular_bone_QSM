function affineMatrix_ijk2RAFplus = get_affineMatrix_ijk2RAFplus(ImDataParams)
% affineMatrix_ijk2RAFplus = get_affineMatrix_ijk2RAFplus(ImDataParams)
    
    voxelSize_mm = ImDataParams.voxelSize_mm;
    pP = ImDataParams.patientPosition;

    T_DICOMijk2RAFplus = get_affineMatrix_DICOMijk2RAFplus(ImDataParams);
    T1 = T_DICOMijk2RAFplus(1:3, 4);
    T_swap = sign(T_DICOMijk2RAFplus(1:3, 1:3));
    additionOffset_mm = -T_swap * voxelSize_mm(:);
    
    T_DICOMijk2RAFplus(1:3, 4) = T1 + additionOffset_mm;
    
    affineMatrix_ijk2RAFplus = T_DICOMijk2RAFplus;
end