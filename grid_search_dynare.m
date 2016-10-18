%% Grid Search with Dynare
% Perhaps useful for finding the best starting point

global COUNT options_
COUNT.total = 0;
COUNT.ss_notfound = 0;
COUNT.ss_neg = 0;
COUNT.ss_found = 0;

%% Define parameter bounds

estimation_init;


%% Set distance between

FF = options_.EST.variables;
LB = options_.EST.LB;
UB = options_.EST.UB;
DIST = options_.EST.DIST;


%% Extended Grid Search - 10 Params
para = {};

for ii = 1:length(FF)
    %list of places to search for each parameter
    % KEY INPUT: how fine should the grid be?
    grid_fineness = 5;
    dist_from_bound = 0.01;
    para{ii} = linspace(LB.(char(FF(ii))) + dist_from_bound*DIST.(char(FF(ii))), UB.(char(FF(ii))) - dist_from_bound*DIST.(char(FF(ii))), grid_fineness);
end

[X1,X2,X3,X4,X5, X6, X7] = ndgrid(para{1}, para{2}, para{3}, para{4}, para{5}, para{6}, para{7});
how_many_grid_points = prod(size(X1))
% fitresult = arrayfun(@(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) fittingfunction10(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10), X1,X2,X3,X4,X5,X6,X7,X8,X9,X10); %run a fitting on every pair fittingfunction(F(J,K), S(J,K))
tic
fitresult = arrayfun(@(p1,p2,p3,p4,p5,p6,p7) fittingfunction10(p1,p2,p3,p4,p5,p6,p7), X1,X2,X3,X4,X5, X6, X7); %run a fitting on every pair fittingfunction(F(J,K), S(J,K))
toc

% I wonder if there are any intersting ways to plot fitresult. 
% Perhaps I start with the 2D case and work my way up. Though each time
% only selecting two dims

% Find the smallest element
[F_min, I] = min(fitresult(:));
F_min
[I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10] = ind2sub(size(fitresult),I)

p1_opt  =  X1(I_row, I_col, I3, I4, I5, I6, I7);
p2_opt  =  X2(I_row, I_col, I3, I4, I5, I6, I7);
p3_opt  =  X3(I_row, I_col, I3, I4, I5, I6, I7);
p4_opt  =  X4(I_row, I_col, I3, I4, I5, I6, I7);
p5_opt  =  X5(I_row, I_col, I3, I4, I5, I6, I7);
p6_opt  =  X6(I_row, I_col, I3, I4, I5, I6, I7);
p7_opt  =  X7(I_row, I_col, I3, I4, I5, I6, I7);
% p8_opt  =  X8(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10);
% p9_opt  =  X9(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10);
% p10_opt = X10(I_row, I_col, I3, I4, I5, I6, I7, I8, I9, I10);

f_opt_real = fittingfunction10(p1_opt, p2_opt, p3_opt, p4_opt, p5_opt, p6_opt, p7_opt)

disp(sprintf('%% Grid Search Results:'))
disp(sprintf([char(FF(1)), ' = %0.10g;'], p1_opt));
disp(sprintf([char(FF(2)), ' = %0.10g;'], p2_opt));
disp(sprintf([char(FF(3)), ' = %0.10g;'], p3_opt));
disp(sprintf([char(FF(4)), ' = %0.10g;'], p4_opt));
disp(sprintf([char(FF(5)), ' = %0.10g;'], p5_opt));
disp(sprintf([char(FF(6)), ' = %0.10g;'], p6_opt));
disp(sprintf([char(FF(7)), ' = %0.10g;'], p7_opt));
% disp(sprintf([char(FF(8)), ' = %0.10g;'], p8_opt));
% % disp(sprintf([char(FF(9)), ' = %0.10g;'], p9_opt));
% disp(sprintf([char(FF(10)), ' = %0.10g;'], p10_opt));

% TODO: what percent of fitresult is 1.0e+10? ie, what percent of parameter
% space gives errors?

% TODO: try with more grid points
% TODO: try with csminwel() being called for each grid point (slow. so only
% do this with very few grid points)

% COUNT

fitresult(fitresult == 10000000000) = NaN;

% surf( X1(:,:,10), X2(:,:,10), fitresult(:,:,10))

% save('grid_search_results_10182016.mat', 'fitresult', 'X1','X2','X3','X4','X5','X6','X7')

% sum(isnan(fitresult(:)))
% ans =
%        22125
% sum(~isnan(fitresult(:)))
% ans =
%        56000

% Timing: A grid of 5^7 (78125) only took 51 minutes to compute

fail_rate = 22125/78125;
success_rate = 56000/78125;

       