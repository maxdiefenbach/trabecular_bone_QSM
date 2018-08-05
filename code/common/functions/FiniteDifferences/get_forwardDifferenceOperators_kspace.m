function [Ei, Ej, Ek] = get_forwardDifferenceOperators_kspace(matrixSize, voxelSize_mm, doIfftshift)
% [Ei, Ej, Ek] = get_forwardDifferenceOperators_kspace(matrixSize, voxelSize_mm, doIfftshift)

    if nargin < 3
        doIfftshift = false;
    end
    if nargin < 2
        voxelSize_mm = [1, 1, 1];
    end

    FOV_mm = matrixSize(:)' .* voxelSize_mm(:)';
    sampleSize = 1 ./ FOV_mm;
    [Ki, Kj, Kk] = get_basisGrid(matrixSize, sampleSize);
    
    Ni = matrixSize(1);
    Nj = matrixSize(2);
    Nk = matrixSize(3);
    
    Ei = (1 - exp(-2j * pi * Ki / Ni)) ./ sampleSize(1);
    Ej = (1 - exp(-2j * pi * Kj / Nj)) ./ sampleSize(2);
    Ek = (1 - exp(-2j * pi * Kk / Nk)) ./ sampleSize(3);
    
    if doIfftshift
        Ei = ifftshift(Ei);
        Ej = ifftshift(Ej);
        Ek = ifftshift(Ek);
    end

end