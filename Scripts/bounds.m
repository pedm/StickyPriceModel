function [ x_bounded ] = bounds( params )
    % csminwel solves an unconstrained minimization problem.
    % this function converts the unbounded params to a bounded guess
    % Inverse function is boundsINV.m 
        
    % IMPORTANT: Bounds are now defined in estimation_init.m
    
    % This file no longer needs to be edited
    
    global options_
    variables = options_.EST.variables;
    LB = options_.EST.LB;
    UB = options_.EST.UB;
    x_bounded = NaN(size(variables));
    
    % For each parameter we are estimating, set the bounds accordingly
    for iii = 1:length(variables)
        lb = LB.(char(variables(iii)));
        ub = UB.(char(variables(iii)));
        x_bounded( iii ) = mintomod_ab(params( iii ), lb, ub);
    end
    

%     x_bounded(1) = mintomod_ab(params(1), 0.1, 0.9); % eta (if eta is large, then ND becomes larger than YD in ss. Thus CD is negative)
%     x_bounded(2) = mintomod_ab(params(2), 0.7, 0.95); % phi
%     x_bounded(3) = mintomod_ab(params(3), 0.01, 0.5); % lambda_ss
%     x_bounded(4) = mintomod_ab(params(4), 0.0001, 0.99); % rhozeta
%     x_bounded(5) = mintomod_ab(params(5),  .05, 10); % sigmazeta If too small, I get "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate."

%     x_bounded(2) = mintomod_ab(params(2), 0.00001, 0.9999); % gamma
%     x_bounded(4) = mintomod_ab(params(4), 0.005, 100); % lambda_bar
%     x_bounded(5) = mintomod_ab(params(5), 0, 100000000); % psi_N
%     x_bounded(7) = mintomod_ab(params(7), 0.0001, 0.99); % rhozeta2
%     % x_bounded(9) = mintomod_ab(params(9),  .1, 10); % zetabar
%     x_bounded(9)= mintomod_ab(params(9), .01, .99); % rho_lambda

end

