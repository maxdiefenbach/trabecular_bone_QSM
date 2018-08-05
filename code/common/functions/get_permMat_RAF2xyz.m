function S = get_permMat_RAF2xyz(PatientPlacement)

    PatientOrientation = PatientPlacement.PatientOrientation;
    PatientPosition = PatientPlacement.PatientPosition;

% set up head first
    T = zeros(3, 3);
    switch PatientOrientation           
      case 'Supine'
        T = [ 0, -1, 0; ...
              1,  0, 0; ...
              0,  0, 1];
      case 'Prone'
        T = [ 0,  1, 0; ...
             -1,  0, 0; ...
              0,  0, 1];
      case 'Left'
        T = [-1,  0, 0; ...
              0, -1, 0; ...
              0,  0, 1];
      case 'Right'
        T = [ 1, 0, 0; ...
              0, 1, 0; ...
              0, 0, 1];
      otherwise
        error(sprintf('Patient position ''%s'' not known.', PatientOrientation))
    end
    
% flip y and z for feet first
    switch PatientPosition
      case 'HeadFirst'
        S = T;
      case 'FeetFirst'
        S = diag([1, -1, -1]) * T;
      otherwise
        error(sprintf('Patient orientation ''%s''not known.', PatientPosition))
    end

end
