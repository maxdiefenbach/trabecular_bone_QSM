function I =  WFI_subject(IDP_filename, save_nii)
    
    if nargin < 2
        save_nii = false;
    end
    if ~isfield(whos, 'I')
        I = ImDataParamsBMRR(IDP_filename);
    end
    assert(strcmp(I.filename, IDP_filename))
    
    outdir = fileparts(I.filename)

    I.multiseedRG
    I.save_WFIparams(fullfile(outdir, [I.fileID '_WFIparams_' I.WFIparams.method '.mat']));

    I.CSS
    
    I.unwrap_fieldmap

    I.WFIparams = rmfield(I.WFIparams, 'paramsMatrixArray')
    fmat = fullfile(outdir, [I.fileID '_WFIparams_' I.WFIparams.method '_unwrap.mat'])
    I.save_WFIparams(fmat);

    if save_nii
        fnii = strrep(I.filename, '_ImDataParams.mat', ['_fieldmap_unwrap.nii.gz']);
        I.save_array2nii(I.WFIparams.fieldmap_Hz_unwrap, fnii);
    end

end
