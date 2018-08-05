function OutParams = closedFormL2(DataParams, Options)

    RDF_ppm = DataParams.RDF_ppm;
    voxelSize_mm = DataParams.voxelSize_mm;
    B0dir = DataParams.B0dir;
    
    lambda = Options.regularizationParameter;
    
    matrixSize = size(RDF_ppm);
    FOV_mm = matrixSize(:)' .* voxelSize_mm(:)';
    DCoffset = 0;
    D = get_dipoleKernel_kspace(matrixSize, FOV_mm, B0dir, DCoffset);
    D = fftshift(D);
    DtD = abs(D).^2;
    
    [k1, k2, k3] = ndgrid(0:matrixSize(1)-1, ...
                          0:matrixSize(2)-1, ...
                          0:matrixSize(3)-1); % already ifftshift-ed
    E1 = 1 - exp(2i .* pi .* k1 / matrixSize(1));
    E2 = 1 - exp(2i .* pi .* k2 / matrixSize(2));
    E3 = 1 - exp(2i .* pi .* k3 / matrixSize(3));
    EtE = abs(E1).^2 + abs(E2).^2 + abs(E3).^2;
    
    chimap_ppm = real(ifftn(D .* fftn(RDF_ppm) ./ (DtD + lambda * EtE + eps)));

    OutParams.chimap_ppm = chimap_ppm;
    OutParams.regularizationParameter = lambda;
    OutParams.method = 'closedFormL2';

end
