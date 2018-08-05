function ImDataParams = phaseDemodulate_field_Hz(ImDataParams, field_Hz)
% ImDataParams = phaseDemodulate_field_Hz(ImDataParams, field_Hz)
% 
% subtracts field_Hz from phase signal

    signal = ImDataParams.signal;
    TE_s = ImDataParams.TE_s;

    for iTE = 1:length(TE_s)
        TE = TE_s(iTE);
        E = exp(-1j * 2 * pi * field_Hz * TE);
        signal(:, :, :, iTE) = signal(:, :, :, iTE) .* E;
    end

    ImDataParams.signal = signal;

end 