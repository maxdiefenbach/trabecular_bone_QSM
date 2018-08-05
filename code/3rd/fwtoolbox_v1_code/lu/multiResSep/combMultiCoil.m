function [sc, coilSense] = srcCoilCombineGEMR1(s, t, df, ksz, stdDev)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine multi-coil source using an algorithm
%% equivalent to GEM R=1
%%
%% s: input signal: nx * ny * nte * ncoils
%% t: echo time in seconds
%% df: 210Hz for 1.5T, 420Hz for 3T
%% ksz: 20
%% stdDev: coil noise standard deviation
%%
%% Huanzhou Yu
%% 060313
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nx = size(s, 1);
ny = size(s, 2);
nte = size(s, 3);
ncoils = size(s, 4);

%stdDev = stdDev/(sum(stdDev))*ncoils;

[mV, mI] = min(abs(round(t*df) - t*df));
ite = mI;
% ite = 2;

filt_kb = zeros(nx, ny);
filt_kb(nx/2+1 - ksz/2: nx/2 + ksz/2, ny/2+1 - ksz/2: ny/2 + ksz/2) = kaiser(ksz, 2)*kaiser(ksz, 2).';
%filt_kb = ones(nx, ny);
for (i1 = 1: ncoils)

  s(:, :, :, i1) = s(:, :, :, i1)/(stdDev(i1)*ncoils);
  S(:, :, i1) = fft2c(s(:, :, ite, i1)).* filt_kb;

  coilSense(:, :, i1) = ifft2c(S(:, :, i1));

  for (i2 = 1: nte)
      s(:, :, i2, i1) = s(:, :, i2, i1).*conj(coilSense(:, :, i1));
  end
end

coilIntensity = sqrt(sum(abs(coilSense).*abs(coilSense), 3));

% dispMR(abs(coilSense(:, :, 1)), 'cs1');
% dispMR(abs(coilSense(:, :, 2)), 'cs2');
% dispMR(abs(coilSense(:, :, 3)), 'cs3');
% dispMR(abs(coilSense(:, :, 4)), 'cs4');
% dispMR(abs(coilSense(:, :,5)), 'cs5');
% dispMR(abs(coilSense(:, :, 6)), 'cs6');
% dispMR(abs(coilSense(:, :, 7)), 'cs7');
% dispMR(abs(coilSense(:, :, 8)), 'cs8');
% dispMR(abs(coilIntensity), 'ci');
%

sc = sum(s, 4);
for (i1= 1: nte)
 sc(:, :, i1) = sc(:, :, i1)./coilIntensity;
end

