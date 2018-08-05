function A = calculate_chemical_shift_encoding_matrix(algoParams,data)

% Author:  Samir Sharma
% Created: November 2011

GAMMA = 42.58; %in MHz/T
for aa = 1:numel(algoParams.species)
  df = GAMMA * data.FieldStrength * algoParams.species(aa).ppm;
  A(:,aa) = exp(1i*2*pi*data.TE(:)*df) * algoParams.species(aa).relAmps(:); % DH*: removed "-" sign from exponential to conform with toolbox
%  A(:,aa) = exp(-1i*2*pi*data.TE(:)*df) * algoParams.species(aa).relAmps(:);
end


