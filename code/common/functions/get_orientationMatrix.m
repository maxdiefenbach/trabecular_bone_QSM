function rotMat_RAFplus2xyz = get_orientationMatrix(patientPosition, patientOrientation)
    
    T = zeros(3, 3);
    switch patientOrientation           % Head First
      case 'Supine'
        T = [ 0, -1, 0;  1,  0, 0; 0, 0, 1];
      case 'Prone'
        T = [ 0,  1, 0; -1,  0, 0; 0, 0, 1];
      case 'Left'
        T = [-1,  0, 0;  0, -1, 0; 0, 0, 1];
      case 'Right'
        T = [ 1,  0, 0;  0,  1, 0; 0, 0, 1];
    else
        error('Patient position not known.')
    end
    
    if strcmp(patientOrientation, 'Head First')
        rotMat_RAFplus2xyz = T;
    elseif strcmp(patientOrientation, 'Feet First')
        rotMat_RAFplus2xyz = diag([1, -1, -1]) * T;
    else
        error('Patient orientation not known.')
    end

end