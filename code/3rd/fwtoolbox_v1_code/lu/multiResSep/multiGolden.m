function [p, rMin, tpsi] = multiGolden(psiRange, s, A, TE, tol, vlevel)
tpsi = []; minResd = [];                
[dpsi, mrd] = golden('rfunc', -psiRange, -1/2*psiRange, 0, ...
    tol,s,A,TE, vlevel);
tpsi = [tpsi dpsi]; minResd = [minResd mrd];
[dpsi, mrd] = golden('rfunc', -1/2*psiRange, 0, 1/2*psiRange, ...
    tol,s,A,TE, vlevel);
tpsi = [tpsi dpsi]; minResd = [minResd mrd];
[dpsi, mrd] = golden('rfunc', 0, 1/2*psiRange, psiRange,  ...
    tol,s,A,TE, vlevel);
tpsi = [tpsi dpsi]; minResd = [minResd mrd];
[rMin, mInd] = min(minResd);

p = tpsi(mInd);
return;
