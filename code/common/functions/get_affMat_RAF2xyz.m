function affMat_RAF2xyz = get_affMat_RAF2xyz(ImageOrientation, longitudinalTableOffcentre_mm)
    
    if nargin < 2
        longitudinalTableOffcentre_mm = 0;
    end

    affMat_RAF2xyz = eye(4);
    
    offcentre_ARF = ImageOrientation.Offcentre
    R_AFR2RAF = get_rotMat_AFR2RAF;
    offcentre_RAF = (R_AFR2RAF * offcentre_ARF(:)).'

    R_RAF2xyz = get_permMat_RAF2xyz(ImageOrientation)
    offcentre_xyz = (R_RAF2xyz * offcentre_RAF(:)).'
    
    offcentre_xyz(3) = offcentre_xyz(3) + longitudinalTableOffcentre_mm

    affMat_RAF2xyz(1:3, 1:3) = get_rotMat_RAF2xyz(ImageOrientation);
    affMat_RAF2xyz(1:3, 4) = offcentre_xyz;

end