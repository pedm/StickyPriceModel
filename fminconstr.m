function [c,ceq] = fminconstr( x )
% wrapper function used in ss solver

    c = [];
    ceq = fbnd(x);

end

% Root Finding with Bounds - Solution Suggested by Matlab:::::::::::::::::::::::::::::::::::::::
% Set Equations and Inequalities as fmincon Constraints
% You can reformulate the problem and use fmincon as follows:
% 
% Give a constant objective function, such as @(x)0, which evaluates to 0 for each x.
% Set the fsolve objective function as the nonlinear equality constraints in fmincon.
% Give any other constraints in the usual fmincon syntax.
% For this example, write a function file for the nonlinear inequality constraint.
% 
% function [c,ceq] = fminconstr(x)
% 
% c = []; % no nonlinear inequality
% ceq = fbnd(x); % the fsolve objective is fmincon constraints
% Save this code as the file fminconstr.m on your MATLAB path.
% 
% Solve the constrained problem.
% 
% lb = [0,0]; % lower bound constraint
% rng default % reproducible initial point
% x0 = 100*randn(2,1);
% x0 = [1.01, 0.5];
% opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
% x = fmincon(@(x)0,x0,[],[],[],[],lb,[],@fminconstr,opts)