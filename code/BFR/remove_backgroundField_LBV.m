function OutParams = remove_backgroundField_LBV(DataParams, Options)
% OutParams = remove_backgroundField_LBV(DataParams, Options)
    
    totalRDF_ppm = DataParams.totalRDF_ppm;
    matrixSize = size(totalRDF_ppm);
    voxelSize_mm = DataParams.voxelSize_mm;
    B0dir = DataParams.B0dir;
    BFRmask = Options.BFRmask;
    %   iFreq - total magnetic field
    %   Mask - ROI
    %   matrix_size - dimension of the 3D image stack
    %   voxel_size - dimensions of the voxels 
    %   tol - iteration stopping criteria on the coarsest grid
    %   depth - number of length scales. The largest length scale is 2^depth * voxel size.
    %   peel - number of boundary layers to be peeled off
    %   N1 - iterations on each depth before the recursive call
    %   N2 - iterations on each depth after the recursive call
    %   N3 - iterations on the finest scale after the FMG is finished.
    %
    %   When using the code, please cite 
    %   Zhou et al. NMR in Biomed 27 (3), 312-319, 2014
    % 
    %   Created by Dong Zhou (zhou.dong@gmail.com) on 2013.06.12
    %   Last modified by Dong Zhou on 2013.06.24
    tol = set_option(Options, 'tol', 0.01);
    depth = set_option(Options, 'depth', -1);
    peel = set_option(Options, 'peel', 0);
    N = set_option(Options, 'nIter', [30, 100, 100]);
    N1 = N(1);
    N2 = N(2);
    N3 = N(3);
    tic
    localRDF_ppm = LBV(totalRDF_ppm, BFRmask, matrixSize, voxelSize_mm, tol, depth, peel, N1, N2, N3);
    elapsedTime_s = toc;

    OutParams.localRDF_ppm = localRDF_ppm;
    OutParams.backgroundRDF_ppm = totalRDF_ppm - localRDF_ppm;
    OutParams.elapsedTime_s = elapsedTime_s;
    OutParams.BFRmethod = 'LBV';
    OutParams = merge_Structs(OutParams, Options);

end