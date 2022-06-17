function sol = Solver(A)
% Parameters
    Beta = 0.5;
    Mu = 0.2;
    [N, ~] = size(A);
    
% Variables
    syms X [1 N];
    X = X.';
    
% Functions definition
    f = simplify(Beta*diag(1-X)*A*X-Mu*X);
    
% Jacobian of f
    J = simplify(jacobian(f, X));
    
x0 = 0.6*ones(N, 1);
[sol,its,err] = Gradient_descent(N, f, J, x0);
sol = round (vpa(sol), 4);
disp('iterations:')
disp(its)
disp('error:')
disp(vpa(err))
disp('zeroes values:')
disp(sol)
end