function imDataParams4toolbox = get_imDataParams4toolbox(ImDataParams)
    
    imDataParams4toolbox = ImDataParams;

    imDataParams4toolbox.images = add_coildim(imDataParams4toolbox.signal);
    imDataParams4toolbox.TE = imDataParams4toolbox.TE_s(:)';
    imDataParams4toolbox.FieldStrength = imDataParams4toolbox.fieldStrength_T;
    imDataParams4toolbox.PrecessionIsClockwise = imDataParams4toolbox.precessionIsClockwise;

    % remove fields
    imDataParams4toolbox = rmfield(imDataParams4toolbox, 'signal');
    imDataParams4toolbox = rmfield(imDataParams4toolbox, 'TE_s');
    imDataParams4toolbox = rmfield(imDataParams4toolbox, 'fieldStrength_T');
    imDataParams4toolbox = rmfield(imDataParams4toolbox, 'precessionIsClockwise');
end


function signal_coils = add_coildim(signal);
    shape = size(signal);
    % nx, ny, nz, ncoils, nTE]
    newShape = [shape(1), shape(2), shape(3), 1, shape(4)];
    signal_coils = reshape(signal, newShape);
end
