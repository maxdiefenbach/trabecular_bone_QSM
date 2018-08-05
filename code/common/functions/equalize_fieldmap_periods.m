function fieldmap_Hz_equalized = equalize_fieldmap_periods(fieldmap_Hz, TE_s)
    
    dTE_s = diff(TE_s(1:2));

    fieldmap_Hz_wrapped = wrap_fieldmap_Hz(fieldmap_Hz, TE_s);
    
    period = round(mean(mean((fieldmap_Hz - fieldmap_Hz_wrapped) * dTE_s, 1), 2));
    
    sz = size(fieldmap_Hz);
    len = 1 / dTE_s;
    fieldmap_Hz_equalized = fieldmap_Hz - repmat(period, [sz(1), sz(2), 1]) * len;

end
