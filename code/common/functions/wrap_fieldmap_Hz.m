function fieldmap_Hz = wrap_fieldmap_Hz(fieldmap_Hz, TE_s)
    
    dTE_s = diff(TE_s(1:2));

    fieldmap_Hz = angle(exp(2j * pi * fieldmap_Hz * dTE_s)) /...
        (2 * pi * dTE_s);

end