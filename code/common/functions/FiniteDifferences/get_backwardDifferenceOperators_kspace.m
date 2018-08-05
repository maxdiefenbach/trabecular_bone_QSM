function [Ei, Ej, Ek] = get_backwardDifferenceOperators_kspace(matrixSize, voxelSize_mm, doIfftshift)
% [Ei, Ej, Ek] = get_backwardDifferenceOperators_kspace(matrixSize, voxelSize_mm, doIfftshift)

    if nargin < 2
        voxelSize_mm = [1, 1, 1];
    end
    
    [Ei, Ej, Ek] = get_forwardDifferenceOperators_kspace(matrixSize, voxelSize_mm, doIfftshift);
    
    Ei = -conj(Ei);
    Ej = -conj(Ej);
    Ek = -conj(Ek);

end