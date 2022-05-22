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
% Jacobian_inverse of f
 Jinv = @(xy) (1/((1+Lam*xy(2))*(1+Lam*xy(1))-(Lam)^2*(1-xy(1))*(1-xy(2))))*[(-1-Lam*xy(1)), (Lam*(xy(1)-1)); 
            (Lam*(xy(2)-1)), (-1-Lam*xy(2))];
% Initial guesses
xy = [0.8; 0.8];
xy = newton2d(f,Jinv,xy);

x = xy(1); y = xy(2);

disp('zeroes values:')
disp(['x = ',num2str(x), '  y = ', num2str(y)])

disp('function values:')
disp(f(xy))
