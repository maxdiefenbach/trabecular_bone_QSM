function rotMat_REC2RAF = get_rotMat_REC2RAF(Angulation_AFR)
    
    R_AFR2RAF = get_rotMat_AFR2RAF;
    Angulation_RAF = (R_AFR2RAF * Angulation_AFR(:)).';
    Z1_Y2_Z3 = rotz(Angulation_RAF(1)) * roty(Angulation_RAF(2)) * rotz(Angulation_RAF(3));
    rotMat_REC2RAF = Z1_Y2_Z3;
end