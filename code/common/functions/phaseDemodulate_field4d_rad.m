function ImDataParams = phaseDemodulate_field4d_rad(ImDataParams, field4d_rad)

    signal = ImDataParams.signal;
    nTE = size(signal, 4);

    for iTE = 1:nTE
        phasor = exp(-1j * field4d_rad(:, :, :, iTE));
        signal(:, :, :, iTE) = signal(:, :, :, iTE) .* phasor;
    end

    ImDataParams.signal = signal;
end