function res = D_t(y)

% adjoint first order finite differences operator
% Mariya Doneva
% May 2009
res = Dx_t(y(:,:,1)) + Dy_t(y(:,:,2));
return;


function res = Dx_t(y)
res = [zeros(size(y,1),1) y(:,1:end-1)] - [y(:,1:end-1) zeros(size(y,1),1)];

function res = Dy_t(y)
res = [zeros(1,size(y,2));y(1:end-1,:)] - [y(1:end-1,:); zeros(1,size(y,2))];


