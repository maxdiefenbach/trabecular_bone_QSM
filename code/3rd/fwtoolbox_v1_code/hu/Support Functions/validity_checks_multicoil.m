fprintf ('\nOriginal Data Matrix is %dx%d, %d slices, %d echoes\n\n', Nx, Ny, Nz, N');

temp=input('Want to reconstruct a particular slice? (0)-No, (>=1)-Yes: ');
if temp>=1  
    temp1 = input('Enter slice range start:end --> ');
    algoParams.sliceofint = temp1;             
else
    algoParams.sliceofint = 1:Nz;               
end
datafull = imDataParams.images;
clear temp;

if mod(algoParams.rg(2),2)==0
    fprintf ('\nPlease select an odd value for the region growing extrapolation kernel size (algoparams.rg(2)');
    fprintf ('\n--> PROGAM ABORTED <--\n');
    break;
end   