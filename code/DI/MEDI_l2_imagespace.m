function QSMparams = MEDI_l2_imagespace(DataParams, Options)
    
    doZeropadding = set_option(Options, 'doZeropadding', false);
    usePreconditioning = set_option(Options, 'usePreconditioning', false);

    % extract input
    voxelSize_mm = DataParams.voxelSize_mm(:)';
    B0dir = DataParams.B0dir;
    psi = DataParams.RDF_ppm;
    dataWeighting = Options.dataWeighting;
    gradWeighting = Options.gradWeighting;
    lambda = Options.regularizationParameter;

    if doZeropadding
        padsize = ceil((size(psi) + 1) / 2);
        psi = padarray(psi, padsize);
        if ~isscalar(dataWeighting)
            dataWeighting = padarray(dataWeighting, padsize);
        end
        if ~isscalar(gradWeighting)
            gradWeighting = padarray(gradWeighting, padsize);
        end
    end

    % fidelity
    kernelSize = size(psi);
    FOV_mm = kernelSize .* voxelSize_mm;
    DCoffset = 0;
    D = get_dipoleKernel_kspace(kernelSize, FOV_mm, B0dir, DCoffset);
    D = ifftshift(D);

    doIfftshift = true;
    [Ei, Ej, Ek] = get_forwardDifferenceOperators_kspace(kernelSize, voxelSize_mm, doIfftshift);
    E2 = conj(Ei) .* Ei + conj(Ej) .* Ej + conj(Ek) .* Ek;
    if usePreconditioning
        Pinv = @(chi) real(ifftn(1 ./ (abs(D).^2 + lambda * E2 + eps) .* fftn(chi)));
    else
        Pinv = @(chi) chi;
    end

    function fid = fidelityTerm(chi)
        fid = real(ifftn(D .* fftn(dataWeighting.^2 .* real(ifftn(D .* fftn(chi))))));
    end
    
    gradWeighting4D = gradWeighting;
    if ndims(gradWeighting) == 3
        gradWeighting4D = repmat(gradWeighting, [1, 1, 1, 3]);
    end

    grad = @(x) fgrad(x);               % use from medi toolbox
    div = @(x) bdiv(x);
    function reg = regularizationTerm(chi)
        reg = div(gradWeighting4D.^2 .* grad(chi));
    end
    
    A = @(chi) Pinv((fidelityTerm(chi) + lambda .* regularizationTerm(chi)));
    b = Pinv(real(ifftn(D .* fftn(dataWeighting.^2 .* psi))));

    Options.CGprecision = set_option(Options, 'CGprecision', 0.01);
    CGoutParams = conjugate_gradients(A, b, Options);
    
    % output
    chimap_ppm = CGoutParams.result;
    residual = CGoutParams.residual;
    if doZeropadding
        chimap_ppm = depad_array3d(chimap_ppm, padsize);
        residual = depad_array3d(residual, padsize);
    end
    QSMparams = rmfield(CGoutParams, {'result'});
    QSMparams.chimap_ppm = chimap_ppm;
    QSMparams.regularizationParameter = lambda;
    QSMparams.residual = residual;
    QSMparams.method = 'MEDIl2TVcg';
end