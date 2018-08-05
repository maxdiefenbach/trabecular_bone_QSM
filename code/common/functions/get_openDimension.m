function openDimension = get_openDimension(ImDataParams)
    
    assert(isfield(ImDataParams, 'orientation'))
    
    orientation = ImDataParams.orientation;
    
    switch orientation
      case 'cor'
        openDimension = 1;
      
      case 'sag'
        openDimension = 2;
      
      case 'tra'
        openDimension = 3;
    end

end