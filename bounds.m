function [ x_bounded ] = bounds( params )
    % csminwel solves an unconstrained minimization problem.
    % this function converts the unbounded params to a bounded guess
    % Inverse function is boundsINV.m 

    % TODO: include gamma, phi, eta, lambda, psi_N, rhozeta, sigmazeta, zetabar

    x_bounded(1) = mintomod_ab(params(1), 0, 1); % eta
    x_bounded(2) = mintomod_ab(params(2), 0, 1); % gamma
    x_bounded(3) = mintomod_ab(params(3), 0, 1); % phi
    x_bounded(4) = mintomod_ab(params(4), 0, 1); % lambda

    x_bounded(5) = mintomod_ab(params(5), 0, 100000000); % psi_N
    x_bounded(6) = mintomod_ab(params(6), 0, 1); % rhozeta
    x_bounded(7) = mintomod_ab(params(7), 0, 100000000); % sigmazeta

end

