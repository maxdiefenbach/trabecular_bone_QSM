function affineMatrix_ijk2xyz = get_affineMatrix_ijk2xyz(ImDataParams)
% affineMatrix_ijk2xyz = get_affineMatrix_ijk2RAFplus(ImDataParams)
    
    voxelSize_mm = ImDataParams.voxelSize_mm;
    iPPmin = ImDataParams.ImagePositionPatientMinMax.imagePositionPatientMin;
    pP = ImDataParams.patientPosition

    T_ijk2RAFplus = get_affineMatrix_ijk2RAFplus(ImDataParams);
    
    T_swap = sign(T_ijk2RAFplus(1:3, 1:3));
    additional_zOffset_mm = T_swap * iPPmin;
    
    R_RAFplus2xyz = get_rotationMatrix_RAFplus2xyz(pP);
    R_RAFplus2xyz(4, 4) = 1
    T_ijk2xyz = R_RAFplus2xyz * T_ijk2RAFplus;
    
    T_RAFplus2xyz = R_RAFplus2xyz;
    T_RAFplus2xyz(4, 4) = 1;
    
    affineMatrix_ijk2xyz = T_ijk2xyz;

end