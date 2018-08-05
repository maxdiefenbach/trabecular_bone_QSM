classdef Conjugate_gradient_Test < matlab.unittest.TestCase
    
    methods (Test)
        
        function test_conjugate_gradient_mat(testCase)
        % example from https://en.wikipedia.org/wiki/Conjugate_gradient_method#Numerical_example
           A = [4, 1; 1, 3]
           b = [1; 2]
           
           xsol = [1/11; 7/11]          % [0.0909, 0.6364]
           
           x = A \ b;
           
           x0 = [2; 1];
           
           tol = 1e-4;
           maxiter = 10;
           verbose = true;
           
           AlgoParams.precision = tol;
           AlgoParams.iterMax = maxiter;
           AlgoParams.verbose = verbose;

           OutParams = conjugate_gradient(A, b, x0, AlgoParams)
           x = OutParams.result
           
           OutParams = conjugate_gradient(A, b, x0)
           x = OutParams.result

           AbsTol = 1e-4;
           testCase.verifyEqual(x, xsol, 'AbsTol', AbsTol);
        end

        
        function test_conjugate_gradient_func(testCase)
        % example from https://en.wikipedia.org/wiki/Conjugate_gradient_method#Numerical_example
           A = [4, 1; 1, 3]
           b = [1; 2]
           
           xsol = [1/11; 7/11]          % [0.0909, 0.6364]
           
           x = A \ b
           
           x0 = [2; 1];
           
           func = @(x) A * x;

           tol = 1e-4;
           maxiter = 10;
           verbose = true;
           
           AlgoParams.precision = tol;
           AlgoParams.iterMax = maxiter;
           AlgoParams.verbose = verbose;

           OutParams = conjugate_gradient(func, b, x0, AlgoParams)
           x = OutParams.result
           
           AbsTol = 1e-4;
           testCase.verifyEqual(x, xsol, 'AbsTol', AbsTol);
        end

    end % methods (Test)

end % EOF
