function rotMat_RAF2xyz = get_rotMat_RAF2xyz(ImageOrientation)
    
    Angulation_AFR = ImageOrientation.Angulation;
    
    R_RAF2xyz = get_permMat_RAF2xyz(ImageOrientation);
    R_AFR2RAF = get_rotMat_AFR2RAF;
    
    Angulation_RAF = (R_AFR2RAF * Angulation_AFR(:)).';
    
    Z1_Y2_Z3 = rotz(Angulation_RAF(1)) * roty(Angulation_RAF(2)) * rotz(Angulation_RAF(3));
    rotMat_RAF2xyz = R_RAF2xyz * Z1_Y2_Z3;
end