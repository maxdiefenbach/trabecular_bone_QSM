function outParams = fw_nsa_map(imDataParams, algoParams, reconParams)



%get water and fat species numbers
waterSpecies = 0; fatSpecies = 0;
for i=1:length(reconParams.species)
    if strcmpi(reconParams.species(1,i).name,'water')
        waterSpecies = i;
    elseif strcmpi(reconParams.species(1,i).name,'fat')
        fatSpecies = i;
    end
end
if ~waterSpecies || ~fatSpecies
    disp('Error - cannot differentiate fat/water species data');
    outParams=[];
    return
end

voxels = size(reconParams.species(1,1).amps);

% DH* Take central slice if multi-slice data
reconParams.species(1).amps = (reconParams.species(1).amps(:,:,ceil(end/2)));
reconParams.species(2).amps = (reconParams.species(2).amps(:,:,ceil(end/2)));
try 
  reconParams.r2starmap = reconParams.r2starmap(:,:,ceil(end/2));
catch
  reconParams.r2starmap = zeros(size(reconParams.species(1).amps(:,:,1)));
end
  
try 
  reconParams.fieldmap = reconParams.fieldmap(:,:,ceil(end/2));
catch
  reconParams.fieldmap = zeros(size(reconParams.species(1).amps(:,:,1)));
end
  


% should first check that voxel #s are consistent across outParams
nsamap.species(waterSpecies).name = 'water';
nsamap.species(waterSpecies).nsa = zeros(voxels);

nsamap.species(fatSpecies).name = 'fat';
nsamap.species(fatSpecies).nsa = zeros(voxels);

%compute NSA at each voxel
gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency)*(imDataParams.FieldStrength)];

% Only single-peak for now, so take the max fat peak for deltaF:
[maxAmp,imaxAmp] = max(abs(algoParams.species(2).relAmps(:)));
dF = gyro*(algoParams.species(2).frequency(imaxAmp)- algoParams.species(1).frequency)*(imDataParams.FieldStrength)
%dF = max(abs(deltaF));



try 
  useR2star = reconParams.include_r2star;
catch
  useR2star = sum(abs(reconParams.r2starmap(:)))>0; % DH*: Use R2* only if recon uses R2*
end

for j=1:voxels(1)
    for k=1:voxels(2)
        waterAmp = abs(reconParams.species(1,waterSpecies).amps(j,k));
        fatAmp = abs(reconParams.species(1,fatSpecies).amps(j,k));
        waterPhase = angle(reconParams.species(1,waterSpecies).amps(j,k));
        fatPhase = angle(reconParams.species(1,fatSpecies).amps(j,k));
        fieldmap = reconParams.fieldmap(j,k);

        
        % DH*: Use R2* only if recon uses R2*
        if useR2star > 0
          r2star = reconParams.r2starmap(j,k);
          %assuming single-peak R2* signal model
          nsa = length(imDataParams.TE)* ...
                (getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap,r2star)./ ...
                 getCRB(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap,r2star));
        else        
          nsa = length(imDataParams.TE)* ...
                (getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap)./ ...
                 getCRB(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap));
          
        end
        
        
        nsaMap.species(waterSpecies).nsa(j,k) = nsa(1);
        nsaMap.species(fatSpecies).nsa(j,k) = nsa(2);
    end
end

%For illustrative and debugginh purposes, we show the map of the NSA for the water image.
% $$$ imagesc(nsaMap.species(waterSpecies).nsa)
% $$$ colormap(gray)
% $$$ colorbar

outParams = nsaMap;

end