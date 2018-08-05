function res = Wavelet_TT(filterType, filterSize, Lx, Ly)

res.adjoint = 0;
res.qmf = MakeONFilter(filterType, filterSize);
res.Lx = Lx;
res.Ly = Ly;
res = class(res,'Wavelet_TT');
