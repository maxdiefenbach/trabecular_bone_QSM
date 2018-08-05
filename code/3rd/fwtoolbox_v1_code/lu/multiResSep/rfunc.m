function [c, WF] = rfunc(psi, S, A, TE)
TE = TE(:);
psiMatInv = diag(exp(-i*2*pi*psi.*TE));

Sc        = psiMatInv*S;
WF        = A\Sc;
c         = norm(Sc-A*WF, 1);