function im2 = dummyCoilCombine3D( im1 )

%% TODO: combine coils as in coilCombine, but in 3D
%% Temporary: use only first coil
[nx,ny,nz,~,N]=size(im1);
im2=reshape(im1(:,:,:,1,:),[nx,ny,nz,N]);

end