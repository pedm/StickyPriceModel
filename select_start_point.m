% Define starting point used in estimation

% Note: to edit the variables used in estimation, edit estimation_init.m

cmd = ['x_start = [' , strjoin(options_.EST.variables, ', '), ']'];
disp(cmd)
eval(cmd)
x_start_unbounded = boundsINV(x_start);