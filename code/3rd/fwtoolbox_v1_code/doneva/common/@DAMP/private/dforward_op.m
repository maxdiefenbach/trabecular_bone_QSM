function dk_data = dforward_op( dx, x, t, f_wf, rel_amp, mask)


dk_data = zeros(size(mask));

water  = x(:,:,1);
fat    = x(:,:,2);
offres = x(:,:,3);

d_water  = dx(:,:,1);
d_fat    = dx(:,:,2);
d_offres = dx(:,:,3);



for j = 1:length(t)
    
    dk_data(:,:,j) = 1i*2*pi*t(j)*(water + fat*sum(rel_amp.*exp(1i*2*pi*f_wf*t(j)))).*exp(1i*2*pi*offres*t(j)).*d_offres + ...
        (d_water + ( d_fat*sum(rel_amp.*exp(1i*2*pi*f_wf*t(j))))).*exp(1i*2*pi*offres*t(j));
    
    dk_data(:,:,j) = mask(:,:,j).*fft2c(dk_data(:,:,j));
    
end

