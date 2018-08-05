function b = getCRB(TE,wm,fm,wp,fp,df,psi,r2s)

    rho = [wm;fm];
    thetas = TE.*df;


    if(nargin < 7) %Dixon model (no field map (or R2*) supplied)
        
        for i=1:length(TE)
            A(2*i-1:2*i,:) = [cos(wp), cos(fp+df*TE(i)*2*pi); ...
                              sin(wp), sin(fp+df*TE(i)*2*pi)];

            A_wp(2*i-1:2*i,:) = [-sin(wp), 0; ...
                                cos(wp), 0];

            A_fp(2*i-1:2*i,:) = [0,-sin(fp+df*TE(i)*2*pi); ...
                                0, cos(fp+df*TE(i)*2*pi)];

        end

        I = zeros(4,4);

        I(1:2,1:2) = A.'*A;
        I(1:2,3) = A.'*A_wp*rho; I(3,1:2) = (A.'*A_wp*rho).';
        I(1:2,4) = A.'*A_fp*rho; I(4,1:2) = (A.'*A_fp*rho).';
        I(3,3) = (A_wp*rho).'*(A_wp*rho);
        I(3,4) = (A_wp*rho).'*(A_fp*rho); 
        I(4,3) = (A_fp*rho).'*(A_wp*rho);
        I(4,4) = (A_fp*rho).'*(A_fp*rho);

        if(cond(I) > 1e16)
            b = Inf*ones(4,1);
        else
            C = inv(I);
            b = [C(1,1);C(2,2);C(3,3);C(4,4)];
        end

        %add given field map, model assuming phases same
        
    elseif(nargin < 8) %IDEAL model (no R2* supplied)
        for i=1:length(TE)
            A(2*i-1:2*i,:) = [cos(wp+psi*TE(i)*2*pi), cos(fp+(df+psi)*TE(i)*2*pi); ...
                              sin(wp+psi*TE(i)*2*pi), sin(fp+(df+psi)*TE(i)*2*pi)];

            A_wp(2*i-1:2*i,:) = [-sin(wp+psi*TE(i)*2*pi), 0; ...
                                cos(wp+psi*TE(i)*2*pi), 0];

            A_fp(2*i-1:2*i,:) = [0,-sin(fp+(df+psi)*TE(i)*2*pi); ...
                                0, cos(fp+(df+psi)*TE(i)*2*pi)];
            
            % added 2pi here
            A_psi(2*i-1:2*i,:) = 2*pi*[-TE(i)*sin(wp+psi*TE(i)*2*pi),-TE(i)*sin(fp+(df+psi)*TE(i)*2*pi); ... % 2 pi missing
                                TE(i)*cos(wp+psi*TE(i)*2*pi), TE(i)*cos(fp+(df+psi)*TE(i)*2*pi)];

        end

        I = zeros(5,5);

        I(1:2,1:2) = A.'*A;
        I(1:2,3) = A.'*A_wp*rho; I(3,1:2) = (A.'*A_wp*rho).';
        I(1:2,4) = A.'*A_fp*rho; I(4,1:2) = (A.'*A_fp*rho).';
        I(1:2,5) = A.'*A_psi*rho; I(5,1:2) = (A.'*A_psi*rho).';
        I(3,3) = (A_wp*rho).'*(A_wp*rho);
        I(3,4) = (A_wp*rho).'*(A_fp*rho); 
        I(4,3) = (A_fp*rho).'*(A_wp*rho);
        I(3,5) = (A_wp*rho).'*(A_psi*rho); 
        I(5,3) = (A_psi*rho).'*(A_wp*rho); 
        I(4,4) = (A_fp*rho).'*(A_fp*rho);
        I(4,5) = (A_fp*rho).'*(A_psi*rho);
        I(5,4) = (A_psi*rho).'*(A_fp*rho);
        I(5,5) = (A_psi*rho).'*(A_psi*rho);
        
        disp('hello')
        F_Pineda = I

        if(cond(I) > 1e18)
            b = Inf*ones(5,1);
        else
            C = inv(I);
            b = [C(1,1);C(2,2);C(3,3);C(4,4);C(5,5)];
        end
        
    else %IDEAL-T2* model (R2* given)
        
        for i=1:length(TE)
            A(2*i-1:2*i,:) = exp(-TE(i)*r2s).* ...
                [cos(wp+psi*TE(i)*2*pi), cos(fp+(df+psi)*TE(i)*2*pi); ...
                 sin(wp+psi*TE(i)*2*pi), sin(fp+(df+psi)*TE(i)*2*pi)];

            A_wp(2*i-1:2*i,:) = exp(-TE(i)*r2s).* ...
                [-sin(wp+psi*TE(i)*2*pi), 0; ...
                 cos(wp+psi*TE(i)*2*pi), 0];

            A_fp(2*i-1:2*i,:) = exp(-TE(i)*r2s).* ...
                [0,-sin(fp+(df+psi)*TE(i)*2*pi); ...
                 0, cos(fp+(df+psi)*TE(i)*2*pi)];
            
            % added 2pi here
            A_psi(2*i-1:2*i,:) = 2 * pi * exp(-TE(i)*r2s).* ... % 2 pi missing
                [-TE(i)*sin(wp+psi*TE(i)*2*pi),-TE(i)*sin(fp+(df+psi)*TE(i)*2*pi); ...
                 TE(i)*cos(wp+psi*TE(i)*2*pi), TE(i)*cos(fp+(df+psi)*TE(i)*2*pi)];

            A_r2s(2*i-1:2*i,:) = exp(-TE(i)*r2s).* ...
                [-TE(i)*cos(wp+psi*TE(i)*2*pi), -TE(i)*cos(fp+(df+psi)*TE(i)*2*pi); ...
                 -TE(i)*sin(wp+psi*TE(i)*2*pi), -TE(i)*sin(fp+(df+psi)*TE(i)*2*pi)];

        end

        I = zeros(6,6);

        I(1:2,1:2) = A.'*A;
        I(1:2,3) = A.'*A_wp*rho; I(3,1:2) = (A.'*A_wp*rho).';
        I(1:2,4) = A.'*A_fp*rho; I(4,1:2) = (A.'*A_fp*rho).';
        I(1:2,5) = A.'*A_psi*rho; I(5,1:2) = (A.'*A_psi*rho).';
        I(1:2,6) = A.'*A_r2s*rho; I(6,1:2) = (A.'*A_r2s*rho).';
        I(3,3) = (A_wp*rho).'*(A_wp*rho);
        I(3,4) = (A_wp*rho).'*(A_fp*rho); I(4,3) = (A_fp*rho).'*(A_wp*rho);
        I(3,5) = (A_wp*rho).'*(A_psi*rho); I(5,3) = (A_psi*rho).'*(A_wp*rho); 
        I(3,6) = (A_wp*rho).'*(A_r2s*rho); I(6,3) = (A_r2s*rho).'*(A_wp*rho); 
        I(4,4) = (A_fp*rho).'*(A_fp*rho);
        I(4,5) = (A_fp*rho).'*(A_psi*rho); I(5,4) = (A_psi*rho).'*(A_fp*rho);
        I(4,6) = (A_fp*rho).'*(A_r2s*rho); I(6,4) = (A_r2s*rho).'*(A_fp*rho);
        I(5,5) = (A_psi*rho).'*(A_psi*rho);
        I(5,6) = (A_psi*rho).'*(A_r2s*rho); I(6,5) = (A_r2s*rho).'*(A_psi*rho);
        I(6,6) = (A_r2s*rho).'*(A_r2s*rho);
        
    if(cond(I) > 1e18)
        b = Inf*ones(6,1);
    else
        C = inv(I);
        b = [C(1,1);C(2,2);C(3,3);C(4,4);C(5,5);C(6,6)];
    end
    
        C = inv(I);                                     
        b = [C(1,1);C(2,2);C(3,3);C(4,4);C(5,5);C(6,6)];

    end
    

end

