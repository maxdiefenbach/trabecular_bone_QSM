function fm_estimate = regiongrowloop (algoParams, ...
Nx,Ny,Nz,xrat,yrat,data,datafull,BW,BWfull,fm_hold)

sup_pix_size = algoParams.rg(1);
extrap_size = algoParams.rg(2);
fm_estimate = zeros(Nx,Ny,Nz);  % Non-zero starting pixel
BWgrow=zeros(Nx,Ny,Nz);         % Binary mask keeping track of region grow

% FOLLOWING Yu's DESCRIPTION
% STEP 4 - Find median and super_pix_size-1 superpixels.

for slicenum = 1:Nz

    fprintf('\n Slice is %d', slicenum);
    med_val_loc = zeros(2,sup_pix_size);
    super = zeros(1,2);
    
    % find and sort ascending all nonzero entries of current slice
    clear temp; temp = BW(:,:,slicenum) .* fm_hold(:,:,slicenum);
    [rr,cc]=find(temp);
    s = zeros(length(rr),1);
    for k = 1:length(rr)
        s(k)=temp(rr(k),cc(k));
    end
    s = sort(s,'ascend'); clear k*;

    %Check, in case s is even, median is averaged
    if mod(length(s),2)==0 
       s = s(1:(length(s)-1));
    end
    smed = median(s);
    
    k=find(s == smed); %find median
    
    med_val = s((k-(sup_pix_size-1)/2):(k+(sup_pix_size-1)/2));
    clear s; clear smed;
    
    % Locate median super pixels coordinates
    kk=0;
    for k = 1:sup_pix_size
        if (length(find (fm_hold(:,:,slicenum)==med_val(k))))==1
            kk = kk+1;
            [med_val_loc(1,kk), med_val_loc(2,kk)] = find(fm_hold(:,:,slicenum) == med_val(k));
        end
    end
    
% STEP 5 - Find COM Super Pixel
    if Nz == 1
        [rr,cc]=cofmass(data(:,:,1).*BW(:,:,1)); %using TE1 weight
    else
        [rr,cc]=cofmass(data(:,:,slicenum,1).*BW(:,:,slicenum,1)); %using TE1 weight
    end
    
% STEP 6 - Find Starting Super Pix closest to COM
    ds = zeros(1,sup_pix_size);
    for k = 1:sup_pix_size
        ds(k)=euclidean([rr,cc],[med_val_loc(1,k),med_val_loc(2,k)]);
    end
    k = find(ds==min(ds));
    if numel(k)>1
	k=k(1);
    end

    super(1)=med_val_loc(1,k);    
    super(2)=med_val_loc(2,k);
    clear temp*; clear ds*;clear k*; clear rr; clear cc;     
        
% STEP 7 - Region Grow
    fprintf ('\n Computed RG seed pixel is %d %d\n', super(1),super(2));

    if Nz == 1        
        temp_mip = abs(data(:,:,1)).*BW;
%        figure;clf;imshow(temp_mip,[]); axis on; %1 slice, so no mip
%        title ('please select your RG seed pixel');
    else
        temp_mip = mip_hh(abs(data(:,:,:,1)).*BW,3); %mipping across z along 1st TE;     
%        figure;clf;imshow(temp_mip,[]); axis on;
%        title ('please select your RG seed pixel');
    end
    %    temp = input (' Input your own RG seed pixel? (0)-No, (>=1)-Yes: ');
    temp = 0;
    if temp >= 1
        [super(2) super(1)]=ginput(1);
        super = round(super);
        fprintf ('\n Your selected RG seed pixel is %d %d', super(1),super(2));
        close all;
    end

    super_val = fm_hold(super(1),super(2),slicenum);        
    super(1) = super(1) * (xrat*2)-1;
    super(2) = super(2) * (yrat*2)-1;
    %med_val_loc(1,:)=med_val_loc(1,:)*(xrat*2)-1;
    %med_val_loc(2,:)=med_val_loc(2,:)*(yrat*2)-1;      
     
% STEP 8 - Define RG Trajectory

    [traj,BWgrow(:,:,slicenum),fm_init,kk] = square_spiral(Nx,Ny,super,super_val,xrat,yrat);
    
    % Run IDEAL on Start Super Pixel Neighborhood.
    clear temp*; [temp1,temp2]=size(fm_init);
    if Nz == 1
        temp_fm = run_ideal(algoParams, temp1,temp2,1,BWfull(kk(1):kk(2),kk(3):kk(4),slicenum), datafull(kk(1):kk(2),kk(3):kk(4),:),fm_init);
    else
        temp_fm = run_ideal(algoParams, temp1,temp2,1,BWfull(kk(1):kk(2),kk(3):kk(4),slicenum), datafull(kk(1):kk(2),kk(3):kk(4),slicenum,:),fm_init);
    end
    fm_estimate(kk(1):kk(2),kk(3):kk(4),slicenum) = temp_fm; % temp contains field map values for super-pixel neighborhood
    clear k*;
  
    fprintf ('\n RG ');
    % NOW REGION GROW
    for kk = 1:length(traj)      
              
        if BWfull(traj(kk,1),traj(kk,2),slicenum)==1
            
            clear temp*
            if mod(kk,floor(length(traj)/100))==0
             fprintf('.' );
            end

            % Define 2D extrapolation kernel and boundary check
            k = extrap_size;
            if ((traj(kk,1)-k)<1)  ; k1x=1;     else k1x = traj(kk,1)-k; end
            if ((traj(kk,1)+k)>Nx) ; k2x=Nx;    else k2x = traj(kk,1)+k; end
            if ((traj(kk,2)-k)<1)  ; k1y=1;     else k1y = traj(kk,2)-k; end
            if ((traj(kk,2)+k)>Ny) ; k2y=Ny;    else k2y = traj(kk,2)+k; end

            % Weighted fm average
            clear fm_new_val;            
            if Nz == 1
                temp2 = BWgrow(k1x:k2x,k1y:k2y,slicenum) .* abs(datafull(k1x:k2x,k1y:k2y,1));
            else
                temp2 = BWgrow(k1x:k2x,k1y:k2y,slicenum) .* abs(datafull(k1x:k2x,k1y:k2y,slicenum,1));
            end
            temp2_sum = sum(sum(temp2));
            temp2 = temp2 / temp2_sum;
            temp3 = BWgrow(k1x:k2x,k1y:k2y,slicenum) .* fm_estimate(k1x:k2x,k1y:k2y,slicenum) .* temp2;            
            fm_new_val = sum(sum(temp3));
            
            if (isnan(fm_new_val))==1                
                fm_new_val=0;
            end
            
            BWgrow(traj(kk,1), traj(kk,2),slicenum)=1;
            
            % ITERATE THIS CURRENT PIXEL AND ASSIGN TO FIELD MAP             
            if Nz == 1
                temp_fm = run_ideal(algoParams, 1,1,1,BWfull(traj(kk,1),traj(kk,2),slicenum),datafull(traj(kk,1),traj(kk,2),:),fm_new_val); 
            else                
                temp_fm = run_ideal(algoParams, 1,1,1,BWfull(traj(kk,1),traj(kk,2),slicenum),datafull(traj(kk,1),traj(kk,2),slicenum,:),fm_new_val);
            end
            fm_estimate(traj(kk,1),traj(kk,2),slicenum) = temp_fm;  
            
        end %if        
    end %kk        
    
    %figure;clf;imshow(fm_estimate(:,:,slicenum),[]);
    
end %slicenum