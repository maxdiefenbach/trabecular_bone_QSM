function res = ifft2c(x)

res = ifftshift(ifft2(fftshift(x)))*sqrt(length(x(:)));
