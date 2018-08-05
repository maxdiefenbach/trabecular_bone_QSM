function water_fat = forward_op_mp1_H(k_data, offres, t, f_wf,rel_amp, mask)
% hermitian forward transform (with fixed offres)
% simulate k-space data for  given fat, water images, offresonance
% 
% 
% x : sizex x sizey x 3 matrix contains water, fat, offres images
% t : echo times 
% mask: sampling mask
% f_wf chemical shift between water and fat

sizex   = size(k_data,1);
sizey   = size(k_data,2);


water_fat = zeros(sizex, sizey,2);

for j = 1:length(t)    
    temp = exp(-1i*2*pi*offres*t(j)).*ifft2c(mask(:,:,j).*k_data(:,:,j)); 
    water_fat(:,:,1) = water_fat(:,:,1) + temp;
    water_fat(:,:,2) = water_fat(:,:,2) + sum(rel_amp.*exp(-1i*2*pi*f_wf*t(j)))*temp;     
end




