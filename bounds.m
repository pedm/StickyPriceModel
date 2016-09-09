function [ x_bounded ] = bounds( params )
    % csminwel solves an unconstrained minimization problem.
    % this function converts the unbounded params to a bounded guess
    % Inverse function is boundsINV.m 

    % TODO: include gamma, phi, eta, lambda, psi_N, rhozeta, sigmazeta, zetabar

    x_bounded(1) = mintomod_ab(params(1), 0, 1); % eta
    x_bounded(2) = params(2); % gamma
    x_bounded(3) = mintomod_ab(params(3), 0, 1); % phi
end

