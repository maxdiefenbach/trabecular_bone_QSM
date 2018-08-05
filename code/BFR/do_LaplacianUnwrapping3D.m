function phaseUnwrapped_rad = do_LaplacianUnwrapping3D(phaseWrapped_rad, voxelSize_mm, doZeropadding)
% phaseUnwrapped_rad = do_LaplacianUnwrapping(phaseWrapped_rad, voxelSize_mm, doZeropadding)
% 
% Schofield and Zhu, 2003 M.A. Schofield and Y. Zhu, Fast phase unwrapping
% algorithm for interferometric applications, Opt. Lett. 28 (2003), pp. 1194?196
    
    if nargin < 3
        doZeropadding = false
    end
    if nargin < 2
        voxelSize_mm = [1, 1, 1]
    end

    matrixSize = size(phaseWrapped_rad);
    if doZeropadding
        padsize = ceil((matrixSize + 1) / 2);
        phaseWrapped_rad = padarray(phaseWrapped_rad, padsize);
        matrixSize = size(phaseWrapped_rad);
    end
    
    doIfftshift = true;
    L = get_LaplaceOperator_kspace(matrixSize, voxelSize_mm, doIfftshift);
    
    Linv = zeros(size(L));
    Linv(L~=0) = 1 ./ L(L~=0);

    A = cos(phaseWrapped_rad) .* ifftn(L .* fftn(sin(phaseWrapped_rad)));
    B = sin(phaseWrapped_rad) .* ifftn(L .* fftn(cos(phaseWrapped_rad)));
    phaseUnwrapped_rad = ifftn(Linv .* fftn(A - B));

    % same as
    % P = exp(1j .* phaseWrapped_rad);
    % phaseUnwrapped_rad = ifftn(Linv .* fftn(imag(conj(P) .* ifftn(L .* fftn(P)))));
     
     if doZeropadding
        phaseUnwrapped_rad = depad_array3d(phaseUnwrapped_rad, padsize);
    end

end