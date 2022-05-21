function xy = newton2d(f,J,xy)
  
    tol = 10^-6;
    err = 1;
    maxits = 100;
    its = 0;
    while (err>tol) && (its<maxits)
         xyold = xy;
         xy = xyold - J(xyold)\f(xyold);
         err = norm(xy-xyold);
         its = its+1;
    end
    
    disp('iterations:')
    disp(its)
    
    disp('error:')
    disp(err)
end