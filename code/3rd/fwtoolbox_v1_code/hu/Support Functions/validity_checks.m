%temp=input('Want to reconstruct a particular slice? (0)-No, (>=1)-Yes: ');
temp = 0;
if temp>=1  
    temp1 = input('Enter slice range start:end --> ');
    algoParams.sliceofint = temp1;             
else
    algoParams.sliceofint = 1:Nz;               
end 

if Nc > 1   
    fprintf ('\n--> MULTI COIL DATA SET <--');
    fprintf ('\nCurrent algorithm will reconstruct only 1 coil at a time');
    temp = input ('\n\nEnter coil # to reconstruct and proceed: ');
    if temp<=Nc
        datafull = imDataParams.images(:,:,:,temp,:);
    else
        fprintf ('--> INVALID ENTRY, PROGRAM ABORTED <-- \n');
        break;
    end
else
    datafull = imDataParams.images;
end              
fprintf ('\nData Matrix is %dx%d, %d slices, %d echoes\n\n', Nx, Ny, Nz, N');
clear temp;

if mod(algoParams.rg(2),2)==0
    fprintf ('\nPlease select an odd value for the region growing extrapolation kernel size (algoparams.rg(2)');
    fprintf ('\n--> PROGAM ABORTED <--\n');
    break;
end   