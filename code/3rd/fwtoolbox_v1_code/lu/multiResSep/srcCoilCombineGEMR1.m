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

[nx, ny, nte, ncoils] = size(s);

%stdDev = stdDev/(sum(stdDev))*ncoils;

[mV, mI] = min(abs(round(t*df) - t*df));
ite = mI;
% ite = 2;

filt_kb = zeros(nx, ny);
filt_kb(nx/2+1 - ksz/2: nx/2 + ksz/2, ny/2+1 - ksz/2: ny/2 + ksz/2) = ...
    kaiser(ksz, 2)*kaiser(ksz, 2).';
%filt_kb = ones(nx, ny);
for (i1 = 1: ncoils)
    
    s(:, :, :, i1) = s(:, :, :, i1)/(stdDev(i1)*ncoils);
    S(:, :, i1) = ft(s(:, :, ite, i1)).* filt_kb; %LPF
    
    coilSense(:, :, i1) = ift(S(:, :, i1));
    
    for (i2 = 1: nte)
        s(:, :, i2, i1) = s(:, :, i2, i1).*conj(coilSense(:, :, i1));
%         s(:, :, i2, i1) = s(:, :, i2, i1).*abs(coilSense(:, :, i1));
  end
end

coilIntensity = sqrt(sum(abs(coilSense).*abs(coilSense), 3));

sc = sum(s, 4);
for (i1= 1: nte)
 sc(:, :, i1) = sc(:, :, i1)./coilIntensity;
end

