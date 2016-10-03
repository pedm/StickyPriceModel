%% Grid Search with Dynare
% Perhaps useful for finding the best starting point

%% Define parameter bounds
% Search over 2^9 grid of parameters

LB.eta = 0;  
LB.gamma = 0.00001;
LB.phi = 0.5;
LB.lambda_bar = 0.005;
LB.psi_N = 0;
LB.rhozeta = 0.0001; 
LB.rhozeta2 = 0.0001;  
LB.sigmazeta = .5;
LB.zetabar = .1;
LB.rho_lambda = .01;

UB.eta = 0.9999;  
UB.gamma = 0.9999; 
UB.phi = 0.99; 
UB.lambda_bar = 10; 
UB.psi_N = 100; 
UB.rhozeta = 1.5; 
UB.rhozeta2 = 0.99; 
UB.sigmazeta = 10;  
UB.zetabar = 10;
UB.rho_lambda = .99;

FF = fields(LB);
 
DIST = [];
DIST.(char(FF(1))) = UB.(char(FF(1))) - LB.(char(FF(1)));
DIST.(char(FF(2))) = UB.(char(FF(2))) - LB.(char(FF(2)));
DIST.(char(FF(3))) = UB.(char(FF(3))) - LB.(char(FF(3)));
DIST.(char(FF(4))) = UB.(char(FF(4))) - LB.(char(FF(4)));
DIST.(char(FF(5))) = UB.(char(FF(5))) - LB.(char(FF(5)));
DIST.(char(FF(6))) = UB.(char(FF(6))) - LB.(char(FF(6)));
DIST.(char(FF(7))) = UB.(char(FF(7))) - LB.(char(FF(7)));
DIST.(char(FF(8))) = UB.(char(FF(8))) - LB.(char(FF(8)));
DIST.(char(FF(9))) = UB.(char(FF(9))) - LB.(char(FF(9)));
DIST.(char(FF(10))) = UB.(char(FF(10))) - LB.(char(FF(10)));



%% Extended Grid Search - 10 Params
para = {};
for ii = 1:length(FF)
    %list of places to search for each parameter
    % KEY INPUT: how fine should the grid be?
    grid_fineness = 2;
    dist_from_bound = 0.1;
    para{ii} = linspace(LB.(char(FF(ii))) + dist_from_bound*DIST.(char(FF(ii))), UB.(char(FF(ii))) - dist_from_bound*DIST.(char(FF(ii))), grid_fineness);
end

% Param 9: search one point only (param 9 is set by ss solver)
ii = 9;
para{ii} = linspace(LB.(char(FF(ii))) + 0.1*DIST.(char(FF(ii))), UB.(char(FF(ii))) - 0.1*DIST.(char(FF(ii))), 1);

[X1,X2,X3,X4,X5,X6,X7,X8,X9,X10] = ndgrid(para{1}, para{2}, para{3}, para{4}, para{5}, para{6}, para{7}, para{8}, para{9}, para{10});
how_many_grid_points = prod(size(X1))
fitresult = arrayfun(@(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) fittingfunction10(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10), X1,X2,X3,X4,X5,X6,X7,X8,X9,X10); %run a fitting on every pair fittingfunction(F(J,K), S(J,K))

% I wonder if there are any intersting ways to plot fitresult. 
% Perhaps I start with the 2D case and work my way up. Though each time
% only selecting two dims

% Find the smallest element
[F_min, I] = min(fitresult(:));
F_min
[I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10] = ind2sub(size(fitresult),I)

p1_opt  =  X1(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10)
p2_opt  =  X2(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10)
p3_opt  =  X3(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10)
p4_opt  =  X4(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10)
p5_opt  =  X5(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10)
p6_opt  =  X6(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10)
p7_opt  =  X7(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10)
p8_opt  =  X8(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10)
p9_opt  =  X9(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10)
p10_opt = X10(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10)

f_opt  = fittingfunction3(p1_opt,p2_opt,p3_opt,p4_opt)
