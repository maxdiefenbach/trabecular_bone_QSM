function grad = compute_central_difference_gradient(array, voxelSize_mm)
% grad = compute_central_difference_gradient(array, voxelSize_mm)
% 
% sets zeros at boundaries

    if nargin < 2
        voxelSize_mm = [1, 1, 1];
    end
    
    switch ndims(array)
      case 3
        grad = compute_central_difference_gradient_3d(array, voxelSize_mm);
      case 2
        grad = compute_central_difference_gradient_2d(array, voxelSize_mm);
      case 1
        grad = compute_central_difference_gradient_1d(array, voxelSize_mm);
    end    
    
end


function grad  = compute_central_difference_gradient_3d(array, voxelSize_mm)
    grad = repmat(zeros(size(array)), [1, 1, 1, 3]);
    
    grad(2:(end-1), :, :, 1) = (array(3:end, :, :) - array(1:(end-2), :, :)) / (2 * voxelSize_mm(1));
    grad(:, 2:(end-1), :, 2) = (array(:, 3:end, :) - array(:, 1:(end-2), :)) / (2 * voxelSize_mm(2));
    grad(:, :, 2:(end-1), 3) = (array(:, :, 3:end) - array(:, :, 1:(end-2))) / (2 * voxelSize_mm(3));
end


function grad  = compute_central_difference_gradient_2d(array, voxelSize_mm)
    grad = repmat(zeros(size(array)), [1, 1, 2]);
    
    grad(2:(end-1), :, 1) = (array(3:end, :) - array(1:(end-2), :)) / (2 * voxelSize_mm(1));
    grad(:, 2:(end-1), 2) = (array(:, 3:end) - array(:, 1:(end-2))) / (2 * voxelSize_mm(2));
end


function grad  = compute_central_difference_gradient_1d(array, voxelSize_mm)
    grad = zeros(size(array));
    
    grad(2:(end-1)) = (array(3:end) - array(1:(end-2))) / (2 * voxelSize_mm(1));
end

