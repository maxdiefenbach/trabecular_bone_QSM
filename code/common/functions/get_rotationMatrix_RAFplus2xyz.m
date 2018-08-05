function rotationMatrix_RAFplus2xyz = get_rotationMatrix_RAFplus2xyz(patientPosition)
% rotationMatrix_RAFplus2xyz = get_rotationMatrix_RAFplus2xyz(patientPosition)
% 
% |--------------------+--------------------+------------+---------------------------------|
% | Patient Position   | Dicom Abbreviation | xyz-System | Affine Transformation Matrix    |
% |--------------------+--------------------+------------+---------------------------------|
% | Head First, Supine | 'HFS'              | PA-RL-FH   | [ 0, -1, 0;  1,  0, 0; 0, 0, 1] |
% | Head First, Prone  | 'HFP'              | AP-LR-FH   | [ 0,  1, 0; -1,  0, 0; 0, 0, 1] |
% | Head First, Left   | 'HFDL'             | LR-PA-FH   | [-1,  0, 0;  0, -1, 0; 0, 0, 1] |
% | Head First, Right  | 'HFDR'             | RL-AP-FH   | [ 1,  0, 0;  0,  1, 0; 0, 0, 1] |
% | Feet First, Supine | 'FFS'              | PA-LR-HF   | [ 0, -1, 0; 1,  0, 0; 0, 0, -1] |
% | Feet First, Prone  | 'FFP'              | AP-RL-HF   | [ 0,  1, 0; 1,  0, 0; 0, 0, -1] |
% | Feet First, Left   | 'FFDL'             | LR-AP-HF   | [-1,  0, 0; 0,  1, 0; 0, 0, -1] |
% | Feet First, Right  | 'FFDR'             | RL-PA-HF   | [ 1,  0, 0; 0, -1, 0; 0, 0, -1] |
% |--------------------+--------------------+------------+---------------------------------|

    rotationMatrix_RAFplus2xyz = nan(3);
    
    switch patientPosition
        
      case 'HFS'
        rotationMatrix_RAFplus2xyz = [ 0, -1, 0;  1,  0, 0; 0, 0, 1];

      case 'HFP' 
        rotationMatrix_RAFplus2xyz = [ 0,  1, 0; -1,  0, 0; 0, 0, 1];

      case 'HFDL'  
        rotationMatrix_RAFplus2xyz = [-1,  0, 0;  0, -1, 0; 0, 0, 1];

      case 'HFDR'
        rotationMatrix_RAFplus2xyz = [ 1,  0, 0;  0,  1, 0; 0, 0, 1];

      case 'FFS'
        rotationMatrix_RAFplus2xyz = [ 0, -1, 0; 1,  0, 0; 0, 0, -1];

      case 'FFP' 
        rotationMatrix_RAFplus2xyz = [ 0,  1, 0; 1,  0, 0; 0, 0, -1];

      case 'FFDL'  
        rotationMatrix_RAFplus2xyz = [-1,  0, 0; 0,  1, 0; 0, 0, -1];

      case 'FFDR' 
        rotationMatrix_RAFplus2xyz = [ 1,  0, 0; 0, -1, 0; 0, 0, -1];
    
    end
    
end % EOF