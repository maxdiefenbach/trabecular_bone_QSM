
function [A C D] = solveA(algoParams)

A = zeros(2*algoParams.N, 2*algoParams.M);
C = zeros(algoParams.N,algoParams.M); 
D = zeros(algoParams.N,algoParams.M);


algoParams.species(1)
algoParams.species(2)



for j = 1:algoParams.M 
   for n = 1:algoParams.N
        C(j,n)= sum(algoParams.species(j).relAmps .* cos(2*pi* algoParams.species(j).frequency * algoParams.te(n)) );
        D(j,n)= sum(algoParams.species(j).relAmps .* sin(2*pi* algoParams.species(j).frequency * algoParams.te(n)) );
   end
   
   A(1:algoParams.N, j*2-1) = C(j,:);
   A(1:algoParams.N, j*2)   = -1*D(j,:);
   A((algoParams.N+1):2*algoParams.N, j*2-1) = D(j,:);
   A((algoParams.N+1):2*algoParams.N, j*2)   = C(j,:);
end
return;