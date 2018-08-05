function dx = dforward_opH(dk_data, x, t, f_wf,rel_amp, mask)



water    = x(:,:,1);
fat      = x(:,:,2);
offres   = x(:,:,3);

d_water  = zeros(size(water));
d_fat    = zeros(size(fat));
d_offres = zeros(size(offres));

dx       = zeros(size(water,1), size(water,2), 3); 
      

    for j = 1:length(t)  
        
        temp     =     ifft2c(mask(:,:,j).*dk_data(:,:,j));        
        
                
        d_water  = d_water   + conj(exp(1i*2*pi*offres*t(j))).*temp;
         
        d_fat    = d_fat     + (sum(rel_amp.*exp(-1i*2*pi*t(j)*f_wf))*conj(exp(1i*2*pi*t(j)*offres))).*temp;     
        
        d_offres = d_offres  + conj(1i*2*pi*t(j)*exp(1i*2*pi*offres*t(j)).*(water +  sum(rel_amp.*exp(1i*2*pi*t(j)*f_wf))*fat)).*temp;  
         
        

    end
    
    dx(:,:,1) = d_water;
    dx(:,:,2) = d_fat;
    dx(:,:,3) = d_offres;
end
