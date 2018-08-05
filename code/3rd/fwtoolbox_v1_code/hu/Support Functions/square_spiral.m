function [trajectory, BW, fm_init, kk] = square_spiral(Nx,Ny,super,super_val,xrat,yrat)
% CREATE SQUARE TRAJECTORY for RG
% =========================================================================
% Nx, Ny - size of bounding matrix
% super - coordinates of starting super pixel, with neighborhood k
% super_val - value of start super pixel
% trajectory - [:,2] sized matrix, holding x and y coordinates of RG traj
% BW - 2D matrix holding starting kernel and neighborhood - set to 1
% fm_init - intialized fm, with correctly allocated 
% kk - bounding coordinates of start super pixel neighborhood

BW = zeros(Nx,Ny); 
xx = super(1);
yy = super(2);
kx = xrat/2; ky = yrat/2;
% kx = 2; ky = 2;

% Boundary check
if ((xx-kx)<1)  ; k1x=1;   else k1x = xx-kx; end
if ((xx+kx)>Nx) ; k2x=Nx;  else k2x = xx+kx; end
if ((yy-ky)<1)  ; k1y=1;   else k1y = yy-ky; end
if ((yy+ky)>Ny) ; k2y=Ny;  else k2y = yy+ky; end
    
BW(k1x:k2x,k1y:k2y)=1;
fm_init = ones(k2x-k1x+1,k2y-k1y+1) * super_val;

% find maximum distance from center to edge
D = max([abs(xx-[1 Nx]) abs(yy-[1 Ny])]); 

count = 1;
trajectory = zeros(numel(BW)-nnz(BW), 2);

k = min(kx,ky);

for kk = (k+1):D
    
    k1x = xx-kk; 
    if (k1x<1)  
        k1x = 1;   
    end
    k2x = xx+kk; 
    if (k2x>Nx) 
        k2x = Nx;  
    end
    k1y = yy-kk; 
    if (k1y<1)  
        k1y = 1;   
    end
    k2y = yy+kk; 
    if (k2y>Ny) 
        k2y = Ny;  
    end

    for r = (k1x+1):k2x
            if (BW(r,k1y)==0)
                trajectory(count,1)=r; trajectory(count,2)=k1y;
                count = count+1;
                BW(r,k1y)=1;
            end
    end

    for c = (k1y+1):k2y
            if (BW(k2x,c)==0)
                trajectory(count,1)=k2x; trajectory(count,2)=c;
                count = count+1;
                BW(k2x,c)=1;
            end
    end

    for r = (k2x-1):-1:k1x
            if (BW(r,k2y)==0)
                trajectory(count,1)=r; trajectory(count,2)=k2y;
                count = count+1;    
                BW(r,k2y)=1;
            end
    end

    for c = (k2y-1):-1:k1y
            if (BW(k1x,c)==0)
                trajectory(count,1)=k1x; trajectory(count,2)=c;
                count = count+1;
                BW(k1x,c)=1;
            end
    end
end   

% Boundary check
if ((xx-kx)<1)  ;k1x=1;   else k1x = xx-kx; end
if ((xx+kx)>Nx) ;k2x=Nx;  else k2x = xx+kx; end
if ((yy-ky)<1)  ;k1y=1;   else k1y = yy-ky; end
if ((yy+ky)>Ny) ;k2y=Ny;  else k2y = yy+ky; end     
BW = zeros(Nx,Ny);
BW(k1x:k2x,k1y:k2y)=1;
kk = [k1x,k2x,k1y,k2y];

%save trajectory.mat trajectory kk;