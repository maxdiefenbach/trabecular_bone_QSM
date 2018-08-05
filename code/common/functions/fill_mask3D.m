function maskFilled = fill_mask3D(mask3D, openDimension)
% maskFilled = fill_mask3D(mask3D)
% maskFilled = fill_mask3D(mask3D, openDimension)

    mask3D = logical(mask3D);

    if nargin < 2
        maskFilled = imfill(mask3D, 'holes');
    
    else
        matrixSize = size(mask3D);
        
        breadSize = matrixSize;
        breadSize(openDimension) = 1;
        bread = ones(breadSize);

        sandwich = cat(openDimension, bread, mask3D, bread);
        
        maskFilled = imfill(sandwich, 'holes');
        
        maskFilled = get_subArray(maskFilled, openDimension, 2:(matrixSize(openDimension)+1));
    end
    
    if islogical(mask3D)
        maskFilled = logical(maskFilled);
    end

end