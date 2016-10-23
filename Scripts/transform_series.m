function [ output_series ] = transform_series( simul )
% Perform same transformations used in PVAR in Stata
% Take logs, take first differences, then apply Helmert's transform

log_series = log( simul );

% replace any cases where R_simul < 0 (prevent Inf and imag numbers)
log_series(simul <= 0) = NaN;

% FD and Helmert Transform
fd_series = first_difference( log_series );
output_series = helmert(fd_series);

end

