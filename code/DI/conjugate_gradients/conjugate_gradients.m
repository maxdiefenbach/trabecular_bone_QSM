function OutParams = conjugate_gradients(A, b, AlgoParams, x0)
% OutParams = conjugate_gradients(A, b, AlgoParams, x0)
%
% iteratively solve A x = b by the conjugate gradient method
% 
% input:
%       - A: matrix A * x = b or function handel A(x) = b
%       - b: array
%       - x0: array, start value for the iteration
%       - AlgoParams.CGprecision:
%           algorithm stops if 
%           norm(A x - b) / norm(b) < CGprecision
%       - AlgoParams.iterMax:
%           maximal number of iterations
%       - AlgoParams.iterRestart:
%           number of iterations before CG restart
% output:
%       - OutParams.result: array, approximative solution x_N to A x = b
%                           after last iteration step
%       - OutParams.totalPrecision:  residual A x_N - b
%                                    after last iteration step N
%       - OutParams.totalPrecision:  norm(A x_N - b) / norm(b) 
%                                    after last iteration step N
%       - OutParams.iterations: number of iterations N
%       
% formulation from the book
% Jonathan Richard Shewchuck: 
% An Introduction to the Conjugate Gradient Method Without the Agonizing Pain
% Appendix B2
    
    AisaFunctionHandle = isa(A, 'function_handle');

    % default
    CGprecision = 1e-3;
    iterMax = 200;
    verbose = 1;
    iterRestart = 50;

    if nargin < 4
        x0 = zeros(size(b));
    end

    % extract AlgoParams
    if nargin > 2
        if isfield(AlgoParams, 'CGprecision')
            CGprecision = AlgoParams.CGprecision;
        end
        if isfield(AlgoParams, 'iterMax')
            iterMax = AlgoParams.iterMax;
        end
        if isfield(AlgoParams, 'iterRestart')
            iterRestart = AlgoParams.iterRestart;
        end
        if isfield(AlgoParams, 'verbose')
            verbose = AlgoParams.verbose;
        end
        if isfield(AlgoParams, 'precision')
            CGprecision = AlgoParams.precision;
        end

    end
    
    matrixSize = size(x0);

    % define abbreviations and assign initial values
    b = b(:);
    x = x0(:);
    
    if AisaFunctionHandle
        Ax = A(reshape(x, matrixSize));
        Ax = Ax(:);
    else
        Ax = A * x;
    end
    r = b - Ax;

    d = r;
    rTr = r' * r;                       % delta_new
    bTb = b' * b;                       % delta_0
    
    totalPrecision = sqrt(rTr / bTb);
    bestPrecision = totalPrecision;
    bestPrecisionIteration = 0;
    iter = 0;
    if verbose
        fprintf('CG: start...\n');
    end
    tic
    while ((iter < iterMax) & (totalPrecision > CGprecision))
        
        % update alpha
        if AisaFunctionHandle
            Ad = A(reshape(d, matrixSize));
            Ad = Ad(:);
        else
            Ad = A * d
        end
        alpha = rTr / (d' * Ad);
        
        % update x
        x = x + alpha * d;
        
        % update r
        if mod(iter, iterRestart) == 0 & iter > 0
            if verbose
                fprintf('CG: %d. iteration, restart CG\n', iter);
            end

            if AisaFunctionHandle
                Ax = A(reshape(x, matrixSize));
                Ax = Ax(:);
            else
                Ax = A * x;
            end
            
            r = b - Ax;
        else
            r = r - alpha * Ad;
        end
        
        % update beta
        rTr_old = rTr;
        rTr = r' * r;
        beta = rTr / rTr_old;
        
        % update d
        d = r + beta * d;
        
        % update rest
        iter = iter + 1;
        totalPrecision = sqrt(rTr / bTb);
        if totalPrecision < bestPrecision
            bestPrecision = totalPrecision;
            bestPrecisionIteration = iter;
        end
        
        if isnan(totalPrecision)
            totalPrecision = 1;
        end

        if verbose
            fprintf('CG: iter = %d, totalPrecision = %f bestPrecision = %f (iter %d)\n', ...
                    iter, totalPrecision, bestPrecision, bestPrecisionIteration);
        end
        
%         if totalPrecision > 1
%             fprintf('CG: abort.\n')
%             break
%         end

    end
    
    time = toc;
    if verbose
        fprintf('CG: ...done\n')
        fprintf('CG: elapsed time: %f seconds \n', time)
        fprintf('CG: #iterations = %d, totalPrecision = %f\n', iter, totalPrecision);
    end
    
    OutParams.result = reshape(x, matrixSize);
    OutParams.residual = reshape(r, matrixSize);
    OutParams.totalPrecision = totalPrecision;
    OutParams.bestPrecision = bestPrecision;
    OutParams.bestPrecisionIteration = bestPrecisionIteration;
    OutParams.iterations = iter;
    OutParams.time_s = time;
end

