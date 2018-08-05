%% Function name: regionGrow
%%
%% Description: Multi-seeded region growing scheme described in:
%%
%% Berglund J, Johansson L, Ahlström H, Kullberg J. Three-point Dixon method enables whole-body 
%% water and fat imaging of obese subjects. Magn Reson Med. 2010, Jun;63(6):1659-1668.
%%
%% Input:
%%   - A: complex linear model matrix, size[3,2]
%%   - bA: complex field map phasor image 1st candidate, size[nx,ny,nz]
%%   - bB: complex field map phasor image 2nd candidate, size[nx,ny,nz]
%%   - c1: threshold on magnitude weight for seed points, size[1]
%%   - c2: threshold on |log(W/F)| for seed points, size[1]
%%   - mw: real-valued magnitude weight image, size[nx,ny,nz]
%%   - voxelSize: voxel dimensions in mm, size[3]
%%
%% Output:
%%   - b: complex field map phasor image, size[nx,ny,nz]
%%
%% Author: Johan Berglund
%% Date created: November 16, 2011
%% Date last modified: November 18, 2011

function b = regionGrow(A,bA,bB,c1,c2,mw,voxelSize)
    res=true(size(bA)); %true=bA, false=bB
    determined=false(size(bA)); %keeps track on determined voxels

    %% Find seed points
    fg=mw>c1; %seed points must have magnitude weight >c1

    d1=(A(1,1)*A(2,2)-A(2,1)*A(1,2))*A(3,2);
    d2=(A(2,1)*A(3,2)-A(3,1)*A(2,2))*A(1,2);
    d3=(A(2,1)*A(3,2)-A(3,1)*A(2,2))*A(1,1);
    d4=(A(1,1)*A(2,2)-A(2,1)*A(1,2))*A(3,1);
    
    %|log(RA)|, RA calculated equivalently to eq. 7
    scoreA=abs(log((d1*bA-d2*bB)./(d3*bB-d4*bA)));
    %|log(RB)|, RB calculated equivalently to eq. 8
    scoreB=abs(log((d1*bB-d2*bA)./(d3*bA-d4*bB)));

    %voxels with mw>c1 and |logR|<c2 are seeds
    determined(fg)=or(scoreA(fg)<c2,scoreB(fg)<c2); 
    %pick solution with smallest |logR| in seeds
    res(determined)=scoreA(determined)<scoreB(determined);

    %% Do region growing
    %create status image; 0=undetermined, 1=bA, 2=bB
    status=zeros(size(bA),'uint8');
    status(and(determined,res))=1;
    status(and(determined,~res))=2;

    %call subroutine RG written in c++
    b_angle = RG(status, angle(bA), angle(bB), mw, voxelSize);
    
    %% Return result as a phasor image
    b = exp(complex(0,b_angle));
end