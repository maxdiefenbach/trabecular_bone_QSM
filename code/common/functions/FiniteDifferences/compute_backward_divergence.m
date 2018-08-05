function div = compute_backward_divergence(grad, voxelSize_mm)
    
    if nargin < 2
        voxelSize_mm = [1, 1, 1];
    end
    
    switch ndims(grad)
      case 4
        div = compute_backward_divergence_3d(grad, voxelSize_mm);
      case 3
        div = compute_backward_divergence_2d(grad, voxelSize_mm);
      case 2
        div = compute_backward_divergence_1d(grad, voxelSize_mm);
    end    
    
end


function div  = compute_backward_divergence_3d(grad, voxelSize_mm)
    div = zeros([size(grad, 1), size(grad, 2), size(grad, 3)]);
    
    div(2:end, :, :) = div(2:end, :, :) - squeeze(grad(2:end, :, :, 1) - grad(1:(end-1), :, :, 1)) / voxelSize_mm(1);
    div(1, :, :) = squeeze(div(1, :, :)) - squeeze(grad(end, :, :, 1) - grad(1, :, :, 1)) / voxelSize_mm(1);
    
    div(:, 2:end, :) = div(:, 2:end, :) - squeeze(grad(:, 2:end, :, 2) - grad(:, 1:(end-1), :, 2)) / voxelSize_mm(2);
    div(:, 1, :) = squeeze(div(:, 1, :)) - squeeze(grad(:, end, :, 2) - grad(:, 1, :, 2)) / voxelSize_mm(2);

    div(:, :, 2:end) = div(:, :, 2:end) - squeeze(grad(:, :, 2:end, 3) - grad(:, :, 1:(end-1), 3)) / voxelSize_mm(3);
    div(:, :, 1) = squeeze(div(:, :, 1)) - squeeze(grad(:, :, end, 3) - grad(:, :, 1, 3)) / voxelSize_mm(3);
end


function div  = compute_backward_divergence_2d(grad, voxelSize_mm)
    div = zeros([size(grad, 1), size(grad, 2)]);
    
    div(2:end, :, :) = div(2:end, :, :) - squeeze(grad(2:end, :, :, 1) - grad(1:(end-1), :, :, 1)) / voxelSize_mm(1);
    div(1, :, :) = squeeze(div(1, :, :)) - squeeze(grad(end, :, :, 1) - grad(1, :, :, 1)) / voxelSize_mm(1);
    
    div(:, 2:end, :) = div(:, 2:end, :) - squeeze(grad(:, 2:end, :, 2) - grad(:, 1:(end-1), :, 2)) / voxelSize_mm(2);
    div(:, 1, :) = squeeze(div(:, 1, :)) - squeeze(grad(:, end, :, 2) - grad(:, 1, :, 2)) / voxelSize_mm(2);
end


function div  = compute_backward_divergence_1d(grad, voxelSize_mm)
    div = zeros(size(grad));
    
    div(2:end) = div(2:end) - (grad(2:end) - grad(1:(end-1))) / voxelSize_mm(1);
    div(1) = div(1) - (grad(end) - grad(1)) / voxelSize_mm(1);
end

