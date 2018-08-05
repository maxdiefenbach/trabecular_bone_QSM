function affineMatrix_DICOMijk2RAFplus = get_affineMatrix_DICOMijk2RAFplus(ImDataParams)
% affineMatrix_DICOMijk2RAFplus = get_affineMatrix_ijk2RAFplus(ImDataParams)
% 
% Notation from:
% http://nipy.org/nibabel/dicom/dicom_orientation.html#dicom-affines-again
% 
% Note that
% affineMatrix_DICOMijk2RAFplus
% transforms from DICOM pixel array coordinate system which starts counting
% with 0 instead of Matlabs array coordinate system counting from 1!

    voxelSize_mm = ImDataParams.voxelSize_mm;
    iOP = ImDataParams.imageOrientationPatient;
    iPPmin = ImDataParams.ImagePositionPatientMinMax.imagePositionPatientMin;
    iPPmax = ImDataParams.ImagePositionPatientMinMax.imagePositionPatientMax;
    
    matrixSize = get_matrixSize(ImDataParams);
    
    N = matrixSize(3);
    T1 = iPPmin;
    TN = iPPmax;
    dr = voxelSize_mm(1);
    dc = voxelSize_mm(2);
    ds = voxelSize_mm(3);
    k = (TN - T1) ./ (N - 1);

    F = [iOP(4:6), iOP(1:3)];
    
    A_multi = eye(4);
    A_multi(1:3, 1:2) = F * diag([dr, dc]);
    A_multi(1:3, 3) = k;
    A_multi(1:3, 4) = T1;
    
    affineMatrix_DICOMijk2RAFplus = A_multi;
end