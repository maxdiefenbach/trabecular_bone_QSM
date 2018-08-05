% function res = rm2d2( maps,  map_init, position, seed, amp_map, non, mode, display)
%
% 2d Region merging algorithm for B0 mapping
%
% This function obtains a B0 field map from a set of possible map values
% Input:
%   - maps: set of possible field map values for each pixel
%   - map_init: initial field map
%   - position: the coordinates of the seed pixel
%   - seed:     the value of the seed pixel
%   - amp_map: map of pixel positions included in the estimation
%   - non: number of neighbours flag
%   - mode:
%           0: selection is based on mean of all valid pixels in the
%           neighbourhood
%           1: selection is based either on either all valid pixels or only
%           pixels selected in the current iteration
%   - display: display flag
%
%
%
% Mariya Doneva, November 2009

function res = rm2d2( maps, map_init, position, seed, amp_map, non, mode, display)


mapsize   = size(maps);
imsize    = mapsize(1:end-1);
no_maps   = mapsize(end);
r_label   = zeros(imsize);

%============== starting point =======================
sp          = position;

if    ( non == 1)
    sp_big    = ones(9,1)*sp;		    % (3^2) Duplicate array, for offset calculations
elseif( non == 2)
    sp_big    = ones(25,1)*sp;		    % (5^2) Duplicate array, for offset calculations
elseif(non == 3)
    sp_big    = ones(49,1)*sp;		    % (7^2) Duplicate array, for offset calculations
elseif(non == 4)
    sp_big    = ones(81,1)*sp;		    % (9^2) Duplicate array, for offset calculations
end

r_label(sp(1),sp(2)) = 1;

% ======= Initialize output variables =======
res               = map_init;          % Allocate B0 map
res(sp(1), sp(2)) = seed;              % set inital map value

% ======= Define coordinate arrays =======
%
xmap = zeros(imsize);		% x-locations
ymap = zeros(imsize); 		% x-locations

% x locations.
m1 = (1:imsize(1))'*ones(1,imsize(2));
xmap(:,:)=m1;

% y locations.
m1 = ones(imsize(1),1)*(1:imsize(2));
ymap(:,:)=m1;


% ======= Compute distances of between all pixels and the seed pixel and
% order pixels with increasing distance
%
p_coords     = [xmap(:) ymap(:)];		                     % Make a list of pixel coordinates.
p_offsets    = p_coords - ones(prod(imsize),1)*sp;           % Make a list of offsets to seed pixel
p_dist       = sqrt(sum(p_offsets.*p_offsets,2));          % Euclidean distances from seed pixel

[pixel_dist,pixel_order] = sort(p_dist); 	                 % Order distances from seed pixel


% Pixel indices starting from upper left corner of the image
x_ind = (1:imsize(1));
y_ind = (1:imsize(2));



% ====== Set neighbourhood region mask ========

if (non == 1)
    % %  Compute offsets to neighbours for a 3x3 neighbourhood
    
    nbr_offsts = [1 1; 2 1; 3 1;
        1 2; 2 2; 3 2;
        1 3; 2 3; 3 3];
    nbr_offsts = 2*ones(9,2)-nbr_offsts;
    
elseif (non == 2)
    % Compute offsets to neighbours for a 5x5 neighbourhood
    
    nbr_offsts = [1 1 ; 2 1 ; 3 1 ; 4 1 ; 5 1 ;
        1 2 ; 2 2 ; 3 2 ; 4 2 ; 5 2 ;
        1 3 ; 2 3 ; 3 3 ; 4 3 ; 5 3 ;
        1 4 ; 2 4 ; 3 4 ; 4 4 ; 5 4 ;
        1 5 ; 2 5 ; 3 5 ; 4 5 ; 5 5 ];
    
    
    nbr_offsts = 3*ones(25,2)-nbr_offsts;
    
elseif(non == 3)
    % Compute offsets to neighbours for a 7x7 neighbourhood
    nbr_offsts = [1 1 ; 2 1 ; 3 1 ; 4 1 ; 5 1 ; 6 1; 7 1;
        1 2 ; 2 2 ; 3 2 ; 4 2 ; 5 2 ; 6 2; 7 2;
        1 3 ; 2 3 ; 3 3 ; 4 3 ; 5 3 ; 6 3; 7 3;
        1 4 ; 2 4 ; 3 4 ; 4 4 ; 5 4 ; 6 4; 7 4;
        1 5 ; 2 5 ; 3 5 ; 4 5 ; 5 5 ; 6 5; 7 5;
        1 6 ; 2 6 ; 3 6 ; 4 6 ; 5 6 ; 6 6; 7 6;
        1 7 ; 2 7 ; 3 7 ; 4 7 ; 5 7 ; 6 7; 7 7];
    nbr_offsts = 4*ones(49,2)-nbr_offsts;
    
elseif(non == 4)
    % Compute offsets to neighbours for a 9x9 neighbourhood
    nbr_offsts = [1 1 ; 2 1 ; 3 1 ; 4 1 ; 5 1 ; 6 1; 7 1; 8 1; 9 1;
        1 2 ; 2 2 ; 3 2 ; 4 2 ; 5 2 ; 6 2; 7 2; 8 2; 9 2;
        1 3 ; 2 3 ; 3 3 ; 4 3 ; 5 3 ; 6 3; 7 3; 8 3; 9 3;
        1 4 ; 2 4 ; 3 4 ; 4 4 ; 5 4 ; 6 4; 7 4; 8 4; 9 4;
        1 5 ; 2 5 ; 3 5 ; 4 5 ; 5 5 ; 6 5; 7 5; 8 5; 9 5;
        1 6 ; 2 6 ; 3 6 ; 4 6 ; 5 6 ; 6 6; 7 6; 8 6; 9 6;
        1 7 ; 2 7 ; 3 7 ; 4 7 ; 5 7 ; 6 7; 7 7; 8 7; 9 7;
        1 8 ; 2 8 ; 3 8 ; 4 8 ; 5 8 ; 6 8; 7 8; 8 8; 9 8;
        1 9 ; 2 9 ; 3 9 ; 4 9 ; 5 9 ; 6 9; 7 9; 8 9; 9 9];
    nbr_offsts = 5*ones(81,2)-nbr_offsts;
    
end



% Start main loop

for k=1:prod(imsize)
    
    % Display map selection
    
    if ((floor(100*(k)/prod(imsize))-floor(100*(k-1)/prod(imsize))>0) || k==1)
        pcdone = floor(100*k/prod(imsize));
        if ((display==1)&& mod(pcdone,5)==0)
            figure(100);
            imshow(real(res(:,:,1).*(res(:,:,1) < 1000000)),[]);title('Field map');
            drawnow;
            
            
            figure(101);
            imshow(r_label,[]);title('Selected region');
            drawnow;
            
        end;
    end;
    
    pixel_num = pixel_order(k);
    
    %include only pixels within the mask
    if (amp_map(x_ind(xmap(pixel_num)) , y_ind(ymap(pixel_num)))==1)
        
        % ===== Extract possible B0 map candidates from maps        
        
        nbr_data   = reshape(maps(x_ind(xmap(pixel_num)) , y_ind(ymap(pixel_num)) ,:),1,no_maps);
        mapdistvec = zeros(1,no_maps);
        
        % ===== Compute neighbors coordinates
        
        if ( non == 1)
            f = 1:9;
        elseif(non ==2)
            f = 1:25;
        elseif(non ==3)
            f = 1:49;
        elseif(non ==4)
            f = 1:81;
        end
                       
        
        cands     = nbr_offsts(f,:);
        ndelts    = cands + ones(length(f),1)*p_offsets(pixel_num,:) ;
        ncoords   = ndelts + sp_big(1:length(f),:);     	% image coordinates of neighbours
        
        ff1 = 0;
        ff2 = 0;
        for p = 1:length(f)
            if((ncoords(p,:)>0) & (ncoords(p,:)<=imsize ))
                if(amp_map(ncoords(p,1),ncoords(p,2))==1)
                    
                    if((r_label(ncoords(p,1),ncoords(p,2)) == 1))
                        ff1(p) = f(p);              % list of already selected neighbours
                    else
                        ff2(p) = f(p);              % list of rest of the neighbours
                    end
                end
            end
        end
        
        
        f1 = ff1(ff1~=0);
        f2 = ff2(ff2~=0);
        
        vtot1 = 0;
        vtot2 = 0;
        
        
        % Compute the average field map in the neighbourhood regions
        
        ncount1 = 0;
        if (~isempty(f1))
            
            
            for p=1:length(f1)
                
                if ( (ncoords(f1(p),:)>0) & (ncoords(f1(p),:)<=imsize )) 
                    vtot1 = vtot1 + res(ncoords(f1(p),1),ncoords(f1(p),2));
                    ncount1 = ncount1 + 1;
                end;
                
            end;
            
            vtot1 = vtot1(:);
        end;
        
        ncount2 = 0;
        if (~isempty(f2))
            
            
            for p=1:length(f2)
                
                if ( (ncoords(f2(p),:)>0) & (ncoords(f2(p),:)<=imsize ))  
                    vtot2 = vtot2 + res(ncoords(f2(p),1),ncoords(f2(p),2));
                    ncount2 = ncount2 + 1;
                end;
                
            end;
            
            vtot2 = vtot2(:);
        end;
        
        
        if (mode ==1)
            if ((ncount1 > pi/4*ncount2)&&(ncount1 + ncount2 > (length(f)-1)/2))
                vtot = vtot1/ncount1;
            else
                vtot = (vtot2 + vtot1)/(ncount1 + ncount2);
            end
        else
            vtot = (vtot2 + vtot1)/(ncount1 + ncount2);
        end
        
       
        
        
        % make sure that at least one neighbour is included
             
        if(vtot~=0)
            % Compute distance from possible field map values to mean
            % neighborhood value
            for mapindex = 1: no_maps
                
                mapdistvec(mapindex) = (nbr_data(mapindex) - vtot)^2;
                
            end;
            
            % Choose the value with minimal distance to the average field map
            [i1 , i2] = min(mapdistvec);
                       
        
        % ===== Store B0 map for the next iteration
        
        res(x_ind(xmap(pixel_num)) , y_ind(ymap(pixel_num))) = nbr_data(i2);
        r_label(x_ind(xmap(pixel_num)) , y_ind(ymap(pixel_num))) = 1;                                  % label pixel as already estimated
        end;
         
    end;
    
end;









