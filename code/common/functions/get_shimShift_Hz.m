function shimShift_Hz = get_shimShift_Hz(ImDataParams)
    
    GYRO = 42.5774806e6;   % gyromagnetic ratio of the proton, [Hz/T]
    LamorFreq_Hz = GYRO * ImDataParams.fieldStrength_T;
    shimShift_Hz = ImDataParams.centerFreq_Hz - LamorFreq_Hz;

end