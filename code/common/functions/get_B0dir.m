function B0dir = get_B0dir(ImDataParams)
    
    voxelSize_mm = ImDataParams.voxelSize_mm;
    iOP = ImDataParams.imageOrientationPatient;
    iPPmin = ImDataParams.ImagePositionPatientMinMax.imagePositionPatientMin;
    iPPmax = ImDataParams.ImagePositionPatientMinMax.imagePositionPatientMax;
    pp = ImDataParams.patientPosition;
    
    matrixSize = get_matrixSize(ImDataParams);
    
    N = matrixSize(3);
    T1 = iPPmin;
    TN = iPPmax;
    dr = voxelSize_mm(1);
    dc = voxelSize_mm(2);
    k = (TN - T1) / (N - 1);

    F = [iOP(4:6), iOP(1:3)];
    Fp = [iOP(4:6) * dr, iOP(1:3) * dc];

    Amulti = eye(4);
    Amulti(1:3, 1:2) = Fp;
    Amulti(1:3, 3) = k;
    Amulti(1:3, 4) = T1;
    
    B0dir = Amulti(1:3, 1:3) \ [0; 0; 1]; % TODO check sign
    B0dir = B0dir / norm(B0dir);

end