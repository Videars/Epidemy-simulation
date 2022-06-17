% Parameters
    Beta = 0.5;
    Mu = 0.2;
    Lam = Beta/Mu;
    N = 2;
% Matrix import
      A_hom = randi(2,N,N) - 1;
      A_hom = A_hom - tril(A_hom,-1) + triu(A_hom,1)';
      A_hom = A_hom - diag(diag(A_hom));
%     A_hom = cell2mat(struct2cell(load('A_hom.mat')));
%     A_sm = cell2mat(struct2cell(load('A_sm.mat')));
%     A_sf = cell2mat(struct2cell(load('A_sf.mat')));

if det(A_hom) ~= 0
% Variables
    disp(det(A_hom));
    syms x1 x2
    X = [x1; x2];
    
% Functions definition
    f_hom = simplify(Beta*diag(1-X)*A_hom*X-Mu*X);
    disp("generazione funzione");
    disp(f_hom);
%     f_sm = Beta*diag(1-X)*A_sm*X-Mu*X;
%     f_sf = Beta*diag(1-X)*A_sf*X-Mu*X;

% Jacobian of f
    J_hom = simplify(jacobian (f_hom, X));
    disp("jacobian definition");
    disp(J_hom);
%     J_sm = jacobian (f_sm, X);
%     J_sf = jacobian (f_sf, X);

% Jacobian inverse of f
    Jinv_hom = simplify(inv(J_hom));
    disp("Inverse jacobian");
    disp(Jinv_hom);
%     Jinv_sm = simplify(inv(J_sm));
%     Jinv_sf = simplify(inv(J_sf));

% Initial guesses
    x0 = 0.6*ones(N, 1);
    sol_hom = vpa(newton(f_hom, Jinv_hom, x0));
%     sol_sm = vpa(newton(f_sm, Jinv_sm, x0));
%     sol_sf = vpa(newton(f_sf, Jinv_sf, x0));

disp('zeroes values for hom:')
disp(sol_hom)
disp('function values for f_hom:')
disp(subs(f_hom, X, sol_hom))
end
% disp('zeroes values for sm:')
% disp(sol_sm)
% disp('function values for f_sm:')
% disp(subs(f_sm, X, sol_sm))
% 
% disp('zeroes values for sf:')
% disp(sol_sf)
% disp('function values for f_sf:')
% disp(subs(f_sf, X, sol_sf))
