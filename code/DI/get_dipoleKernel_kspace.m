function dipoleKernel_kspace = get_dipoleKernel_kspace(matrixSize, FOV_mm, B0dir, DCoffset)
% dipoleKernel_kspace = get_dipoleKernel_kspace(matrixSize, FOV_mm, B0dir, DCoffset)
% 
% construct dipole kernel in k-space
    
    if nargin < 4
        DCoffset = 0;
    end

    [I, J, K] = get_basisGrid(matrixSize);
    
    samplingInterval_ijk = 1 ./ FOV_mm;
    Ki = I * samplingInterval_ijk(1);
    Kj = J * samplingInterval_ijk(2);
    Kk = K * samplingInterval_ijk(3);
    
    Kz = Ki * B0dir(1) + Kj * B0dir(2) + Kk * B0dir(3);
    K2 = Ki.^2 + Kj.^2 + Kk.^2;
    
    dipoleKernel_kspace = 1/3 - Kz.^2 ./ K2;
    
    dipoleKernel_kspace(isnan(dipoleKernel_kspace)) = DCoffset;
    
end
