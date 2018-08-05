load trajectory.mat;

trajectory_map = zeros(Nx,Ny);
trajectory_map(kk(1):kk(2),kk(3):kk(4))=1;
[m,n]=size(trajectory);

figure(1);clf;imshow(trajectory_map,[]);
title ('region growing trajectory');
hold on;

for k = 1:m
    trajectory_map(trajectory(k,1),trajectory(k,2))=1;
    imshow(trajectory_map,[]); drawnow;
end    

% TOTAL TRAJECTORY COVERAGE
% trajectory_map(trajectory(:,1),trajectory(:,2))=1;
% figure(2);clf;imshow(trajectory_map,[]);