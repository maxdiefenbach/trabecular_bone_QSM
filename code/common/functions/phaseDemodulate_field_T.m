function ImDataParams = phaseDemodulate_field_T(ImDataParams, field_T)
    
    GYRO = 42.5774806e6;   % gyromagnetic ratio of the proton, [Hz/T]
    field_Hz = GYRO * field_T;
    ImDataParams = phaseDemodulate_field_Hz(ImDataParams, field_Hz);
    ImDataParams.isDemodulated = true;

end
