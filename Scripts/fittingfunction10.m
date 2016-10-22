function [ f ] = fittingfunction10( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)

global options_

% Pass these parameters to the distance function

% Currently only works for 3 variables
x_start = [p1, p2, p3, p4, p5, p6, p7];
x_start_unbounded = boundsINV(x_start);
check_bounds;
f = distance_fcn(x_start_unbounded);
    
end


