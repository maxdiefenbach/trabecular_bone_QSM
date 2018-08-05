function cx = cls4(A, b)

% MATLAB code for complex least squares approximation with constrained phase.
% Finds min[norm(A*x*exp(i*phi)-b)] over real vector x and real scalar phi.
%
% Reference:
% Markovsky, I. (2011)
% On the complex least squares problem with constrained phase.
% SIAM Journal on Matrix Analysis and Applications, 32 (3). pp. 987-992. 

C = [real(A) real(b) -imag(b); 
     imag(A) imag(b)  real(b)];
n = size(A, 2); D = diag([zeros(1, n), 1, 1]);  
R   = triu(qr(C, 0)); 
[u, s, v] = svd(R((n + 1):(n + 2), (n + 1):end));
phi = angle(v(1, 2) - i * v(2, 2)); 
x   = R(1:n, 1:n) \ (R(1:n, (n + 1):end) * [v(1, 2); v(2, 2)]);
cx  = x * exp(i * phi);
