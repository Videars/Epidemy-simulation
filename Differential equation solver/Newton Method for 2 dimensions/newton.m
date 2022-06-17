function Xi = newton(f,Jinv,Xi)
  
    tol = 10^-6;
    err = 1;
    maxits = 100;
    its = 0;
    syms x1 x2
    X = [x1; x2];
    
    while (err>tol) && (its<maxits)
         Xold = Xi;
         Xi = Xold - vpa((subs(Jinv, X, Xold)))*vpa((subs(f, X, Xold)));
         disp(Xi);
         err = norm(Xi-Xold);
         its = its+1;
         disp(its);
    end
    
    disp('iterations:')
    disp(its)
    
    disp('error:')
    disp(vpa(err))
end
