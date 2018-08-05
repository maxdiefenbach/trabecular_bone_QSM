function I = make_nii_subject(IDP_filename, WFImethodCell, BFRmethodCell, QSMmethodCell)

    if nargin < 4
        QSMmethodCell = {'closedFormL2', 'MEDIl2TVcg', 'MEDIl1TVnesta'};
    end
    if nargin < 3
        BFRmethodCell = {'LBV'};
    end
    if nargin < 2
        WFImethodCell = {'CSS_unwrap'};
    end

    
    if ~isfield(whos, 'I')
        I = ImDataParamsBMRR(IDP_filename);
    end
    assert(strcmp(I.filename, IDP_filename))

    MIP = I.get_echoMIP;
    fnii = strrep(IDP_filename, 'ImDataParams.mat', 'echoMIP.nii.gz');
    I.save_array2nii(MIP, fnii);


    for WFImethodC = WFImethodCell
        WFImethod = WFImethodC{1}
        I.load_WFIparams(strrep(I.filename, 'ImDataParams', ...
                                ['WFIparams_' WFImethod]))


        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['R2s_' WFImethod '.nii.gz']);
        I.save_array2nii(I.WFIparams.R2s_Hz, fnii);

        I.set_fatFraction_percent;
        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['PDFF_' WFImethod '.nii.gz']);
        I.save_array2nii(I.WFIparams.fatFraction_percent, fnii);
        
        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['fieldmap_' WFImethod '.nii.gz']);
        I.save_array2nii(I.WFIparams.fieldmap_Hz, fnii);

        fnii = strrep(IDP_filename, 'ImDataParams.mat', ['fieldmap_unwrap_' WFImethod '.nii.gz']);
        I.save_array2nii(I.WFIparams.fieldmap_Hz_unwrap, fnii);


        for BFRmethodC = BFRmethodCell
            BFRmethod = BFRmethodC{1}
            I.load_BFRparams(strrep(I.filename, 'ImDataParams', ...
                                    ['BFRparams_' BFRmethod '_' WFImethod]))
            mask3D = binary_erode(I.BFRparams.BFRmask, 1);

            fnii = strrep(IDP_filename, 'ImDataParams.mat', ...
                          ['localRDF_' BFRmethod '_' WFImethod '.nii.gz']);
            I.save_array2nii(I.BFRparams.localRDF_ppm, fnii);


            for QSMmethodC = QSMmethodCell
                QSMmethod = QSMmethodC{1}
                I.load_QSMparams(strrep(I.filename, 'ImDataParams', ...
                                        ['QSMparams_' QSMmethod '_' BFRmethod '_' WFImethod]))

                fnii = strrep(IDP_filename, 'ImDataParams.mat', ...
                              [QSMmethod '_' BFRmethod '_' WFImethod '.nii.gz']);
                I.save_array2nii(mask3D .* I.QSMparams.chimap_ppm, fnii);
                
            end
        end
    end


end

