function R = get_rotMat_dicom2RAI(angulationRAI_deg)
    
    angulationRAI_rad = (pi / 180) .* angulationRAI_deg;
    ang_RL = angulationRAI_rad(1);
    ang_AP = angulationRAI_rad(2);
    ang_IS = angulationRAI_rad(3);

    R_RL = [ 1, 0, 0; ...
             0, cos(ang_RL), -sin(ang_RL);...
             0, sin(ang_RL), cos(ang_RL)];
    
    R_AP = [cos(ang_AP), 0, sin(ang_AP); ...
            0, 1, 0; ...
            -sin(ang_AP), 0, cos(ang_AP)];
    
    R_IS = [cos(ang_IS), -sin(ang_IS), 0; ...
            sin(ang_IS), cos(ang_IS), 0; ...
            0, 0, 1];
    
    R = R_RL * R_AP * R_IS;
    
end