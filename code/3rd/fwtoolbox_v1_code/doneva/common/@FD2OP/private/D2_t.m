function res = D2_t(y)

% adjoint second second order finite differences transform 
% Mariya Doneva
% May 2009
res = D2x_t(y(:,:,1)) + D2y_t(y(:,:,2));

function res = D2x_t(y)
res = - [y(:,1:end-2) zeros(size(y,1),2)] + 2*[zeros(size(y,1),1)  y(:,1:end-2) zeros(size(y,1),1)]  - [zeros(size(y,1),2) y(:,1:end-2) ];   

function res = D2y_t(y)
res = - [y(1:end-2,:);zeros(2,size(y,2))] + 2*[zeros(1,size(y,2)); y(1:end-2,:);zeros(1,size(y,2))]  - [zeros(2,size(y,2)); y(1:end -2,:)];
