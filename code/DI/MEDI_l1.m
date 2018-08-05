function QSMparams = MEDI_l1(DataParams, Options)
    
    Options
    
    % extract input
    psi = DataParams.RDF_ppm;
    voxelSize_mm = DataParams.voxelSize_mm;
    B0dir = DataParams.B0dir;
    
    dataWeighting = Options.dataWeighting;
    gradWeighting = Options.gradWeighting;
    lambda = Options.regularizationParameter;
    mu = Options.smoothingParameter;
    
    % dipole kernel
    kernelSize = size(psi);
    FOV_mm = kernelSize(:)' .* voxelSize_mm(:)';
    DCoffset = 0;
    D = get_dipoleKernel_kspace(kernelSize, FOV_mm, B0dir, DCoffset);
    D = ifftshift(D);

    % prepare nesta
    A = @(x) real(dataWeighting .* ifftn(D .* fftn(x)));
    At = @(y) real(ifftn(D .* fftn(dataWeighting .* y)));
    b = psi;
    op = real(ifftn(D.^2));
    Lipschitz = set_option(Options, 'LipschitzConstant', sqrt(sum(abs(op(:)).^2)));

    gradWeighting4D = gradWeighting;
    if ndims(gradWeighting) == 3
        gradWeighting4D = repmat(gradWeighting, [1, 1, 1, 3]);
    end
    
    opts.gradWeighting4D = gradWeighting4D;
    opts.TolVar = set_option(Options, 'precision', 1e-2);
    opts.Verbose = set_option(Options, 'verbose', 25);
    opts.maxcont = set_option(Options, 'maxcont', 5);
    opts.maxiter = set_option(Options, 'maxiter', 10);
    opts.stopTest = set_option(Options, 'stopTest', 1);
    opts.MaxIntIter = set_option(Options, 'MaxIntIter', 5);
    opts.TypeMin = 'tv';
    opts.voxelSize_mm = voxelSize_mm(:)';
    if isfield(Options, 'chimapInit_ppm')
        opts.xplug = Options.chimapInit_ppm;
        disp('set chimap init');
    end

    fprintf('NESTA: lambda = %f, mu = %f, Lipschitz = %f\n', lambda, mu, Lipschitz)
    fprintf('NESTA: start\n')
    tic;
    [x_nesta, niter, resid, err, optsOut] = NESTA_UP(A, At, b, lambda, Lipschitz, mu, opts);
    elapsedTime_s = toc;

   % output
    QSMparams.chimap_ppm = x_nesta;
    QSMparams.iterations = niter;
    QSMparams.residual = resid;
    QSMparams.error = err;
    QSMparams.elapsedTime_s = elapsedTime_s;
    QSMparams.regularizationParameter = lambda;
    QSMparams.smoothingParameter = mu;
    QSMparams.dataWeighting = dataWeighting;
    QSMparams.method = 'MEDIl1TVnesta';
    QSMparams = merge_Structs(QSMparams, optsOut);
end