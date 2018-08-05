function I = BFR_subject(IDP_filename, WFImethod, save_nii)

    if nargin < 3
        save_nii = false;
    end        
    if nargin < 2
        WFImethod = 'CSS_unwrap';
    end
    
    if ~isfield(whos, 'I')
        I = ImDataParamsBMRR(IDP_filename);
        I.load_WFIparams(strrep(I.filename, 'ImDataParams', ['WFIparams_' WFImethod]))
    end
    assert(strcmp(I.filename, IDP_filename))


    BFRmask = binary_erode(fill_mask3D(I.get_tissueMask(2)), 5);
    
    % LBV
    I.set_BFRparams

    I.BFRparams.BFRmask = BFRmask;
    I.remove_backgroundField_LBV
    fmat = strrep(I.filename, '_ImDataParams.mat', ['_BFRparams_LBV_' WFImethod '.mat']);
    I.save_BFRparams(fmat);

    if save_nii
        fnii = strrep(I.filename, '_ImDataParams.mat', ['_localRDF_LBV_' WFImethod '.nii.gz']);
        I.save_array2nii(I.BFRparams.localRDF_ppm, fnii);
    end
        
% % PDF
% I.BFRparams = []
% I.set_BFRparams

% I.BFRparams.BFRmask = BFRmask;
% I.remove_backgroundField_PDF
% fmat = strrep(I.filename, '_ImDataParams.mat', ['_BFRparams_PDF_' WFImethod '.mat']);
% I.save_BFRparams(fmat);
    
% if save_nii
%     fnii = strrep(I.filename, '_ImDataParams.mat', ['_localRDF_PDF_' WFImethod '.nii.gz']);
%     I.save_array2nii(I.BFRparams.localRDF_ppm, fnii);
% end

end
