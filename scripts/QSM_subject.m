function I = QSM_subject(IDP_filename, WFImethod, BFRmethod, save_nii)

    if nargin < 2
        WFImethod = 'CSS_unwrap';
    end
    if nargin < 3
        BFRmethod = 'LBV';
    end
    if nargin < 4
        save_nii = true;
    end

    if ~isfield(whos, 'I')
        I = ImDataParamsBMRR(IDP_filename);
        I.load_WFIparams(strrep(I.filename, 'ImDataParams', ['WFIparams_' WFImethod]))
        I.load_BFRparams(strrep(I.filename, 'ImDataParams', ['BFRparams_' BFRmethod '_' WFImethod]))
    end
    assert(strcmp(I.filename, IDP_filename))
    
    MIP = I.get_echoMIP;
    I.set_outOfPhaseImg;
    I.set_fatFraction_percent;
    signalMask = threshold_array(MIP, 15);
    mask3D = binary_erode(I.BFRparams.BFRmask, 1);
    I.QSMparams.mask = mask3D;
    if save_nii
        fnii = strrep(IDP_filename, 'ImDataParams.mat', 'echoMIP.nii.gz');
        I.save_array2nii(MIP, fnii);
        
        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['R2s_' WFImethod '.nii.gz']);
        I.save_array2nii(I.WFIparams.R2s_Hz, fnii);

        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['PDFF_' WFImethod '.nii.gz']);
        I.save_array2nii(I.WFIparams.fatFraction_percent, fnii);
        
        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['fieldmap_' WFImethod '.nii.gz']);
        I.save_array2nii(I.WFIparams.fieldmap_Hz, fnii);

        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['fieldmap_unwrap_' WFImethod '.nii.gz']);
        I.save_array2nii(I.WFIparams.fieldmap_Hz_unwrap, fnii);

        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['localRDF_' BFRmethod '_' WFImethod '.nii.gz']);
        I.save_array2nii(I.BFRparams.localRDF_ppm, fnii);
    end

    I.set_QSMparams
    
    %% closedFormL2
    Options = struct();
    Options.regularizationParameter = 0.08;
    I.QSMparams.AlgoParams = Options;
    I.closedFormL2
    fmat = strrep(IDP_filename, 'ImDataParams', ['QSMparams_closedFormL2_' BFRmethod '_' WFImethod]);
    I.save_QSMparams(fmat)
    if save_nii
        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['chi_closedFormL2_' BFRmethod '_' WFImethod '.nii.gz']);
        I.save_array2nii(mask3D .* I.QSMparams.chimap_ppm, fnii);
    end
    
    % dataWeighting
    dataWeighting = scale_array2interval(MIP, [0, 1]);
    
    % gradWeighting
    G = fgrad(MIP);% 4D
    Gsum = sum(abs(G), 4);
    edgeMask4D = 0.1 * max(G(:)) < G;
    gradWeighting4D = ~edgeMask4D;

    %% MEDI l2
    Options = struct();
    Options.regularizationParameter = 0.008;
    Options.dataWeighting = dataWeighting;
    Options.gradWeighting = gradWeighting4D;
    Options.usePreconditioning = 1;
    I.QSMparams.AlgoParams = Options;
    I.MEDI_l2
    fmat = strrep(IDP_filename, 'ImDataParams', ['QSMparams_MEDIl2TVcg_' BFRmethod '_' WFImethod]);
    I.save_QSMparams(fmat)
    if save_nii
        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['chi_MEDIl2TVcg_' BFRmethod '_' WFImethod '.nii.gz']);
        I.save_array2nii(mask3D .* I.QSMparams.chimap_ppm, fnii);
    end

    %% NESTA MEDI l1
    Options = struct();
    Options.regularizationParameter = 0.008;
    Options.smoothingParameter = 1e-6;       
    Options.dataWeighting = dataWeighting;
    Options.gradWeighting = gradWeighting4D;
    I.QSMparams.AlgoParams = Options;
    I.MEDI_l1
    fmat = strrep(IDP_filename, 'ImDataParams', ['QSMparams_MEDIl1TVnesta_' BFRmethod '_' WFImethod]);
    I.save_QSMparams(fmat)
    if save_nii
        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['chi_MEDIl1TVnesta_' BFRmethod '_' WFImethod '.nii.gz']);
        I.save_array2nii(mask3D .* I.QSMparams.chimap_ppm, fnii)
    end

end
