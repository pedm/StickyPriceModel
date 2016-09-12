function [ x_unbounded ] = boundsINV( params )
    % csminwel solves an unconstrained minimization problem.
    % this function converts the bounded params to an unbounded guess
    % Inverse function is bounds.m 

    % TODO: include gamma, phi, eta, lambda, psi_N, rhozeta, sigmazeta, zetabar

    x_unbounded(1) = modtomin_ab(params(1), 0, 1); % eta
    x_unbounded(2) = modtomin_ab(params(2), 0, 1); % gamma
    x_unbounded(3) = modtomin_ab(params(3), 0, 1); % phi
    x_unbounded(4) = modtomin_ab(params(4), 0, 1); % lambda
    
    x_unbounded(5) = modtomin_ab(params(5), 0, 100000000); % psi_N
    x_unbounded(6) = modtomin_ab(params(6), 0, 1); % rhozeta
    x_unbounded(7) = modtomin_ab(params(7), 0, 100000000); % sigmazeta

end

