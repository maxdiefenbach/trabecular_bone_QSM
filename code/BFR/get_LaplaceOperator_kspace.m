function L = get_LaplaceOperator_kspace(matrixSize, voxelSize_mm, doIfftshift)
% L = get_LaplaceOperator_kspace(matrixSize, voxelSize_mm, doIfftshift)
% 
% contruct Laplace kernel in image space and apply fftn

    if nargin < 3
        doIfftshift = false;
    end
   
    kernelSize_ispace = [3, 3, 3];
    kernel_ispace = zeros(kernelSize_ispace);
    h = voxelSize_mm;
    kernel_ispace(:, :, 1) = [0, 0       , 0; ...
                              0, 1/h(3)^2, 0; ...
                              0, 0       , 0];
    kernel_ispace(:, :, 2) = [0       , 1/h(1)^2                       , 0; ...
                              1/h(2)^2, -2/h(1)^2 - 2/h(2)^2 - 2/h(3)^2, 1/h(2)^2; ...
                              0       , 1/h(1)^2                       , 0];
    kernel_ispace(:, :, 3) = [0, 0       , 0; ...
                              0, 1/h(3)^2, 0; ...
                              0, 0       , 0];
    
    kernel_ispaceBig = zeros(matrixSize);
    center = ceil((matrixSize + 1) ./ 2);
    kernelSizeHalf_ispace = ceil((kernelSize_ispace + 1) ./ 2);
    kernel_ispaceBig((center(1)-kernelSizeHalf_ispace(1)+1):(center(1)-kernelSizeHalf_ispace(1)+kernelSize_ispace(1)), ...
                     (center(2)-kernelSizeHalf_ispace(2)+1):(center(2)-kernelSizeHalf_ispace(2)+kernelSize_ispace(2)), ...
                     (center(3)-kernelSizeHalf_ispace(3)+1):(center(3)-kernelSizeHalf_ispace(3)+kernelSize_ispace(3))) = kernel_ispace;

    L = fftn(ifftshift(kernel_ispaceBig));
    if ~doIfftshift
        L = fftshift(L);
    end

end