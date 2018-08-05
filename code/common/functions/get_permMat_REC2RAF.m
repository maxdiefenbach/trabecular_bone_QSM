function permMat_REC2RAF = get_permMat_REC2RAF(ImageOrientation)
    
    orientation = ImageOrientation.Orientation;
    
    switch orientation
      case 'COR'
        T = [ 0,  0, -1;  1,  0,  0;  0,  1,  0];
      case 'SAG'
        T = [ 0,  0, -1;  0,  1,  0; -1,  0,  0];
      case 'TRA'
        T = [ 0,  0,  1;  1,  0,  0;  0,  0,  1];
      otherwise
        error(sprintf('Orientation ''%s'' not known.', orientation));
    end
    
    permMat_REC2RAF = T;
end