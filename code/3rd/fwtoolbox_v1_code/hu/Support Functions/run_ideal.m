function fm_hold = run_ideal(algoParams,Nx,Ny,Nz,BW,data,fm_in)
% constants
maxiter = algoParams.maxiter;
N = algoParams.N;
M = algoParams.M;
te = algoParams.te;
A = algoParams.A;
C = algoParams.C;
D = algoParams.D;

if nargin == 6
    fm_init = zeros(Nx,Ny,Nz);
else
    fm_init = fm_in;
end
fm_hold = zeros(Nx,Ny,Nz);

for slicenum = 1:Nz
    [rr,cc]=find(BW(:,:,slicenum));
    fm = zeros(Nx,Ny);
     
    Ims = zeros(Nx,Ny,N);
    if Nz ==1 
        Ims=data;
    else
        Ims(:,:,1:N)=data(:,:,slicenum,:);
    end
    for nn = 1:N     
             Ims(:,:,nn)=Ims(:,:,nn).* exp(-1i*2*pi*te(nn)  .* fm_init(:,:,slicenum));             
    end

    for k = 1:length(rr)                       
         S_hat = [squeeze(real(Ims(rr(k),cc(k),:))).' squeeze(imag(Ims(rr(k),cc(k),:))).'].';
         p_hat = algoParams.Ainv * S_hat;
         B = calcg(p_hat,A,C,D,te,N,M);         
         S_dhat = zeros(N*2,1);
         for n = 1:N
              temp_sum = 0;
              for j = 1:M
                  temp_sum = temp_sum + ( p_hat(j*2-1)*C(j,n) - p_hat(j*2)*D(j,n) );
              end
              S_dhat(n)=S_hat(n)-temp_sum; % Eq. [B3]
              temp_sum = 0;
              for j = 1:M
                  temp_sum = temp_sum + ( p_hat(j*2-1)*D(j,n) + p_hat(j*2)*C(j,n) );
              end
              S_dhat(n+N) = S_hat(n+N)-temp_sum; % Eq. [B4]
         end
         y = (transpose(B)*B) \ transpose(B) * S_dhat;
        
         del_fm = y(1);
         counter = 1; 

         while (abs(del_fm) > algoParams.fm_epsilon)

             fm(rr(k),cc(k)) = fm(rr(k),cc(k)) + del_fm;

             tempS=reshape(Ims(rr(k),cc(k),:),[N 1]).* exp(-1i*2*pi*fm(rr(k),cc(k))*te');
  
             S_hat = [real(tempS).' imag(tempS).'].';
             p_hat = algoParams.Ainv * S_hat;
             B = calcg(p_hat,A,C,D,te,N,M);             
             S_dhat = zeros(N*2,1);
             for n = 1:N
                  temp_sum = 0;
                  for j = 1:M
                      temp_sum = temp_sum + ( p_hat(j*2-1)*C(j,n) - p_hat(j*2)*D(j,n) );
                  end
                  S_dhat(n)=S_hat(n)-temp_sum; % Eq. [B3]
                  temp_sum = 0;
                  for j = 1:M
                      temp_sum = temp_sum + ( p_hat(j*2-1)*D(j,n) + p_hat(j*2)*C(j,n) );
                  end
                  S_dhat(n+N) = S_hat(n+N)-temp_sum; % Eq. [B4]
             end
             y = (transpose(B)*B) \ transpose(B) * S_dhat;
             
             del_fm = y(1);
             counter = counter+1;
                        
             if counter > maxiter
                del_fm=0;
             end;
        
         end % while loop
         
         if counter <= maxiter
             fm(rr(k),cc(k))=fm(rr(k),cc(k))+y(1);
         else             
             fm(rr(k),cc(k)) = 0;
         end
    end %k
    
    fm_hold(:,:,slicenum) = fm+fm_init(:,:,slicenum); 
  end %loop slicenum
  
