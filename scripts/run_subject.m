function run_subject(IDP_filename)

    if ~isfield(whos, 'I')
        I = ImDataParamsBMRR(IDP_filename);
    end
    assert(strcmp(I.filename, IDP_filename))

    I = WFI_subject(IDP_filename)
    I = BFR_subject(IDP_filename)
    I = QSM_subject(IDP_filename)

end
