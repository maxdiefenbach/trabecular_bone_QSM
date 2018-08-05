function [Ei, Ej, Ek] = get_centralDifferenceOperators_kspace(matrixSize, voxelSize_mm)
% [Ei, Ej, Ek] = get_centralDifferenceOperators_kspace(matrixSize, voxelSize_mm)

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
    
    Ei = (exp(2j * pi * Ki / Ni) - exp(-2j * pi * Ki / Ni)) ./ (2 .* sampleSize(1));
    Ej = (exp(2j * pi * Kj / Nj) - exp(-2j * pi * Kj / Nj)) ./ (2 .* sampleSize(2));
    Ek = (exp(2j * pi * Kk / Nk) - exp(-2j * pi * Kk / Nk)) ./ (2 .* sampleSize(3));

    if doIfftshift
        Ei = ifftshift(Ei);
        Ej = ifftshift(Ej);
        Ek = ifftshift(Ek);
    end

end