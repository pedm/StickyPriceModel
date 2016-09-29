function [ x_bounded ] = bounds( params )
    % csminwel solves an unconstrained minimization problem.
    % this function converts the unbounded params to a bounded guess
    % Inverse function is boundsINV.m 

    % TODO: include gamma, phi, eta, lambda, psi_N, rhozeta, sigmazeta, zetabar

    x_bounded(1) = mintomod_ab(params(1), 0, 0.9999); % eta (if eta is large, then ND becomes larger than YD in ss. Thus CD is negative)
    x_bounded(2) = mintomod_ab(params(2), 0.00001, 0.9999); % gamma
    x_bounded(3) = mintomod_ab(params(3), 0.5, 0.99); % phi
    x_bounded(4) = mintomod_ab(params(4), 0.005, 100); % lambda_bar
    x_bounded(5) = mintomod_ab(params(5), 0, 100000000); % psi_N
    x_bounded(6) = mintomod_ab(params(6), 0.0001, 1.5); % rhozeta
    x_bounded(7) = mintomod_ab(params(7), 0.0001, 0.99); % rhozeta2
    x_bounded(8) = mintomod_ab(params(8),  .5, 10); % sigmazeta If too small, I get "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate."
    x_bounded(9) = mintomod_ab(params(9),  .1, 10); % zetabar

    % sigmazeta 
    % If too small, I get 
    % "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate."
end

