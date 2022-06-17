function [Xi,its,err] = Gradient_descent(N, F, jacF, Xi)
  
    tol = 10^-6;
    err = 1;
    maxits = 100;
    its = 0;
    syms X [1 N];
    X = X.';
    jacF_t = jacF.';
    Xold_old = zeros(N,1);
    while (err>tol) && (its<maxits)
         Xold = Xi;
         Fn_Fn1 = (((vpa(subs(jacF_t, X, Xold)))*(vpa(subs(F, X, Xold))))-((vpa(subs(jacF_t, X, Xold_old)))*(vpa(subs(F, X, Xold_old)))));
         gamma_i = vpa(abs((Xold-Xold_old).'*Fn_Fn1)/(norm(Fn_Fn1)^2));
         Xi = Xold - gamma_i * (vpa(subs(jacF_t, X, Xold)))*(vpa(subs(F, X, Xold)));
         err = norm(Xi-Xold);
         Xold_old = Xold;
         its = its+1;
         disp(its);
    end
end