function res = fft2c(x)

res = fftshift(fft2(ifftshift(x)))/sqrt(length(x(:)));

