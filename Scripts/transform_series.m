function [ output_series ] = transform_series( simul )
% Perform same transformations used in PVAR in Stata
% Take logs, take first differences, then apply Helmert's transform


% replace any cases where simul < 0 (prevent Inf and imag numbers in log_series)
simul(simul <= 0) = NaN;

log_series = log( simul );

% FD and Helmert Transform
fd_series = first_difference( log_series );
output_series = helmert(fd_series);

end

