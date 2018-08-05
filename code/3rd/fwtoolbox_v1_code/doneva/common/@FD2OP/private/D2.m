function res = D2(x)

% second order finite difference transform
% Mariya Doneva
% May 2009

 D2x =  2*[x(:,2:end -1)  zeros(size(x,1),2)]  - [x(:,1:end-2)  zeros(size(x,1),2)] - [x(:,3:end)  zeros(size(x,1),2)];   
 D2y =  2*[x(2:end -1,:); zeros(2,size(x,2))]  - [x(1:end-2,:); zeros(2,size(x,2))] - [x(3:end,:) ;zeros(2,size(x,2))]; 

 res = cat(3,D2x,D2y);



