function [I, J, K] = get_basisGrid(matrixSize, voxelSize, center)
    

    [I, J, K] = ndgrid((-matrixSize(1)/2):(matrixSize(1)/2-1), ...
                       (-matrixSize(2)/2):(matrixSize(2)/2-1), ...
                       (-matrixSize(3)/2):(matrixSize(3)/2-1));
    
    % same as 
    % [J, I, K] = meshgrid((-matrixSize(1)/2):(matrixSize(1)/2-1), ...
    %                      (-matrixSize(2)/2):(matrixSize(2)/2-1), ...
    %                      (-matrixSize(3)/2):(matrixSize(3)/2-1));

    if nargin >= 2
        I = I * voxelSize(1);
        J = J * voxelSize(2);
        K = K * voxelSize(3);
    end
    
    if nargin == 3
        I = I - center(1);
        J = J - center(2);
        K = K - center(3);
    end

end