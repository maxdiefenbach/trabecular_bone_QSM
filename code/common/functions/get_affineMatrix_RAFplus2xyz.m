function affineMatrix_RAFplus2xyz = get_affineMatrix_RAFplus2xyz(patientOrientation, patientPosition, xyzOffcentres)
    
    T = get_rotationMatrix_RAFplus2xyz(patientOrientation, patientPosition);
    T(4, 4) = 1;
    T(1:3, 4) = xyzOffcentres;
    
    affineMatrix_RAFplus2xyz = T;

end