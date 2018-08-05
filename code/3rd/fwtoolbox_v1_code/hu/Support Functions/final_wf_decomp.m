function outputParams = final_wf_decomp(algoParams,Nx,Ny,Nz,filtsize,datafull,BWfull,fm)

Ims_final = zeros(Nx,Ny,Nz,algoParams.M); 
fm_filtered = zeros(Nx,Ny,Nz);
for slicenum = 1:Nz      
    fm_filtered(:,:,slicenum) = medfilt2(fm(:,:,slicenum),filtsize,'symmetric');
end

for slicenum = 1:Nz
    Ims = zeros(Nx,Ny,algoParams.N);
    if Nz == 1
        Ims=datafull;
    else
        Ims=datafull(:,:,slicenum,:);
    end
    for nn = 1:algoParams.N
        
%         if algoParams.Precession==-1
%             Ims(:,:,nn)=Ims(:,:,nn).* exp(-1i*2*pi*algoParams.te(nn).* fm_filtered(:,:,slicenum));
%         else
%             Ims(:,:,nn)=Ims(:,:,nn).* exp(1i*2*pi*algoParams.te(nn).* fm_filtered(:,:,slicenum));
%         end            

             Ims(:,:,nn)=Ims(:,:,nn).* exp(-1i*2*pi*algoParams.te(nn).* fm_filtered(:,:,slicenum));
    end
    
    [rr,cc]=find(BWfull(:,:,slicenum));
    for k = 1:length(rr)   
         S_hat = zeros(algoParams.N*2,1);
         S_hat(1:algoParams.N) = real(Ims(rr(k),cc(k),:)); 
         S_hat((algoParams.N+1):(2*algoParams.N)) = imag(Ims(rr(k),cc(k),:));         
         p_hat = algoParams.Ainv * S_hat;
         for j = 1:algoParams.M
             Ims_final(rr(k),cc(k),slicenum,j) = p_hat((j*2-1)) + 1i*p_hat((j*2));
         end
    end
end

outputParams.water = Ims_final(:,:,:,1);
outputParams.fat = Ims_final(:,:,:,2);
outputParams.fieldmap = fm;

outputParams.species(1).amps = Ims_final(:,:,:,1);
outputParams.species(2).amps = Ims_final(:,:,:,2);