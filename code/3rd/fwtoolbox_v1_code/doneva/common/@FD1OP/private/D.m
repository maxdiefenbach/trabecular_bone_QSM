function res = D(x)

% first order finite differences operator
% Mariya Doneva
% May 2009

Dx  = [x(:,2:end) zeros(size(x,1),1)] - [x(:,1:end-1)  zeros(size(x,1),1)];
Dy  = [x(2:end,:);zeros(1,size(x,2))] - [x(1:end-1,:); zeros(1,size(x,2))];

res = cat(3,Dx,Dy);
