function [newNx,newNy,Nz,data,datafull,BW,BWfull,xrat,yrat] ...
    = getlowresdata (datafull,algoParams)

[Nx Ny , ~, Nc N]=size(datafull);

%% GET DESIRED SLICES
datafull = datafull(:,:,algoParams.sliceofint,:,:);
Nz = numel(algoParams.sliceofint);
%% COIL COMBINE
if Nc > 1   
    fprintf ('\n--> MULTI COIL DATA SET <--');
    fprintf ('\nCombining Data via Diego Hernando coilCombine Function');
    temp = zeros(Nx,Ny,Nz,1,N);
    
%         for sl = 1:Nz
%             temp(:,:,sl,1,:)=coilCombine(datafull(:,:,sl,:,:));
%         end
    temp = coilCombine3D(datafull);   
    clear datafull; datafull = temp; clear temp;
    fprintf ('\nCoil Combine Done');
    Nc = 1;
end         

%% DOWNSAMPLE
newNx = algoParams.downsize(1);
newNy = algoParams.downsize(2);

xrat = round(Nx/newNx/2);
yrat = round(Ny/newNy/2);

if mod(xrat,4)~=0
    xrat = 4;
    newNx = Nx/2/xrat;
end
if mod(yrat,4)~=0
    yrat = 4;
    newNy = Ny/2/yrat;
end

% full ksp data
tempfull = zeros(Nx,Ny,Nz,Nc,N);

% IFT
for nn = 1:N
    tempfull(:,:,:,:,nn)=fftshift(ifft2(fftshift(datafull(:,:,:,:,nn))));               
end 

lowx_resrange = ((Nx-newNx)/2+1) : ((Nx-newNx)/2+newNx);
lowy_resrange = ((Ny-newNy)/2+1) : ((Ny-newNy)/2+newNy);

data = zeros(newNx,newNy,Nz,Nc,N);
temp = tempfull(lowx_resrange,lowy_resrange,:,:,:);

% FT
for nn = 1:N
    data(:,:,:,:,nn) = fftshift(fft2(fftshift(temp(:,:,:,:,nn))));   
end

% remove coil dimension
datafull = squeeze(datafull);
data = squeeze(data);

fprintf ('\n');
%temp1 = input ('Want to draw mask ROI? (0)-No, (>=1)-Yes: ');
temp1 = 0;
if temp1 >= 1
%    subplot(122);imshow(mip_hh(abs(datafull(:,:,:,1)),3),[]); axis on;
    title ('Draw a bounding ROI by clicking mouse anchor points.  Double click to finish draw.  Double click again to continue');
    BWfull = roipoly; 
    BWfull = repmat(BWfull,[1,1,Nz]);
    fprintf ('\nROI Done\n'); 
else
    BWfull = ones(Nx,Ny,Nz); 
end

BW = ones(newNx,newNy,Nz);
% BW = zeros(newNx,newNy,Nz);
% for slicenum = 1:Nz
%     BW(:,:,slicenum) = imresize(BWfull(:,:,slicenum),[newNx,newNy]);
% end
fprintf ('\nData LowRes Matrix is %dx%d, %d slices, %d echoes\n\n', newNx, newNy, Nz, N');