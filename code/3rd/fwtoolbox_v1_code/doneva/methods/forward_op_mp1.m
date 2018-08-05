function k_data = forward_op_mp1(x, t, f_wf, rel_amp, mask)
% forward transform for multipeak fat model
% simulate k-space data for  given fat, water images, offresonance
% 
% 
% x : sizex x sizey x 3 matrix contains water, fat, offres images
% t : echo times 
% mask: sampling mask
% f_wf chemical shift between water and fat
% rel_amp relative amplitude of fat peaks


sizex   = size(x,1);
sizey   = size(x,2);

k_data  = zeros(sizex, sizey,length(t));


water   = x(:,:,1);
fat     = x(:,:,2);
offres  = x(:,:,3);

for j = 1:length(t)     
    k_data(:,:,j) = mask(:,:,j).*fft2c(exp(1i*2*pi*offres*t(j)).*(water + fat*sum(rel_amp.*exp(1i*2*pi*f_wf*t(j)))));   
end

