function B = calcg(p,A,C,D,te,N,M)
g = zeros(2*N,1);
for n = 1:N
     temp_sum = 0;
     for j = 1:M
         temp_sum = temp_sum + ( -1*p(j*2-1)*D(j,n) - p(j*2)*C(j,n) );
     end
     temp_sum = temp_sum * (2*pi*te(n));
     g(n) = temp_sum;         
end

for n = (N+1):(2*N)
     temp_sum = 0;
     for j = 1:M
         temp_sum = temp_sum + ( p(j*2-1)*C(j,n-N) - p(j*2)*D(j,n-N) );
     end
     temp_sum = temp_sum * (2*pi*te(n-N));
     g(n) = temp_sum;              
end
B = [g,A];