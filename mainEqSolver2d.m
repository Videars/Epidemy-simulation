%Parameters
 Beta=0.5;
 Mu=0.2;
 Lam=Beta/Mu;

% Functions 
 f = @(xy) [-xy(1) + Lam*(1-xy(1))*xy(2);
            -xy(2) + Lam*(1-xy(2))*xy(1)];   
% Jacobian of f
 J = @(xy) [1 - Lam*xy(2), Lam*(1-xy(1)); 
            Lam*(1-xy(2)), 1 - Lam*xy(1)];
        
% Initial guesses
xy = [0.5; 0.5];
xy = newton2d(f,J,xy);

x = xy(1); y = xy(2);

disp('zeroes values:')
disp(['x = ',num2str(x), '  y = ', num2str(y)])

disp('function values:' )
disp(f(xy))