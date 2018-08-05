
%% Function: createSynthetic_kSpace
%%
%% Description: create synthetic chemical shift-encoded dataset, with arbitrary species and echo times. 
%%
%% Some features:
%%    - k-space
%%    - Single-R2*
%%    - Accepts multi-peak species
%%    - Accepts multiple species
%%    - Accepts a field map in Hz
%%    - single TE for the k-space is assumed 
%%
%% Input arguments:
%%   - kDataParams0.fieldStrength: (in Tesla)
%%   - kDataParams0.kSpaceTimes: (in seconds)
%%   - kDataParams0.kSpaceLocations: locations for each k-space value, array of size[Nsamples,3]
%%
%%   - algoParams.species(ii).name = name of species ii (string)
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0] 
%%      - algoParams.species(1).relAmps = [1]   
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%%
%%   - trueParams.species(ii).amps: true water/fat images, size [nx,ny,nz,ncoils] 
%%   - trueParams.r2starmap: R2* map (in s^{-1}, size [nx,ny,nz])
%%   - trueParams.fieldmap: field map (in Hz, size [nx,ny,nz])
%%
%%
%% Output:
%%   - kDataParams.kSpaceValues: acquired k-space values, array of size[Nsamples,1]
%%   - kDataParams.FOV: (in Tesla)
%%
%%
%% Author: Mariya Doneva
%% Date created: 
%% Date last modified: 
%%


function kDataParams = createSynthetic_kSpace( kDataParams0, algoParams, trueParams )

gyro = 42.58;
kDataParams = kDataParams0;
t = kDataParams0.kSpaceTimes;
[sx,sy,sz,C] = size(trueParams.species(1).amps);

TE = unique(t);     % This is quick fix to avoid computing the k-space one point at a time in case of repeating TE, not recommended for general trajectories
N = length(TE);  

Nsamples = length(t);
try
  fieldStrength = kDataParams0.fieldStrength;
catch
  fieldStrength = 1.5;
end  
  
try 
  r2starmap = trueParams.r2starmap;
catch 
  r2starmap = zeros(sx,sy);
end
  
try 
  fieldmap = trueParams.fieldmap;
catch 
  fieldmap = zeros(sx,sy);
end
   
 kDataParams.kSpaceValues = zeros(Nsamples,1);
 
 kSpaceData = zeros(sx,sy,sz,C,N); 
 
 for kt = 1:N
     
     image = zeros(sx,sy,sz);
     
     for ks=1:length(trueParams.species)
         
         amps = trueParams.species(ks).amps;
         
         try
             freqs = gyro*fieldStrength*algoParams.species(ks).frequency;
         catch
             freqs = 0;
         end
         
         try
             relAmps = algoParams.species(ks).relAmps;
         catch
             relAmps = ones(size(freqs))/length(freqs);
         end
         
                          
         s = sum(relAmps(:).*exp(1i*2*pi*freqs(:)*TE(kt)));
         fieldAndR2starEffect = exp(-TE(kt)*r2starmap + 1i*2*pi*TE(kt)*fieldmap);     
         image = image + amps.*s.*fieldAndR2starEffect;
         
     end
     % compute k-space data
     kSpaceData = fft2c(image);
     locations_idx = find(kDataParams.kSpaceTimes==TE(kt));
    
     idx = sub2ind([sx,sy],kDataParams.kSpaceLocations(locations_idx,1)+ sx/2 +1, kDataParams.kSpaceLocations(locations_idx,2)+ sy/2 +1 );
     temp = kSpaceData(idx);
      
     kDataParams.kSpaceValues(locations_idx) = temp;
           
    
 end
   
 
 kDataParams.FOV = [sx;sy;1];                      % 2x1 double 


