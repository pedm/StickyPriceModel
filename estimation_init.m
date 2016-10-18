%% Initialize Estimation

% List the parameters to estimate
options_.EST = [];
options_.EST.irf_length = 31; % 11 or 31. note: you also have to change this in stoch_simul()
options_.EST.variables = {'eta', 'phi', 'lambda_ss', 'psi_N', 'rhozeta', 'sigmazeta', 'gamma', 'rho_lambda'}; 
% , 'psi_N', 'rhozeta', 'sigmazeta'
% , 'rho_lambda', 'psi_N', 'gamma', 'L_ss'   gamma ?

% to list the parameter values:
% eval(['[' , strjoin(options_.EST.variables, ', '), ']'])

%% Define parameter bounds

options_.EST.LB.eta = 0.1;  
options_.EST.LB.gamma = 0.00001;
options_.EST.LB.phi = 0.7;
options_.EST.LB.psi_N = 0;
options_.EST.LB.rhozeta = 0.0001; 
options_.EST.LB.sigmazeta = .05;
options_.EST.LB.rho_lambda = .01;
options_.EST.LB.lambda_ss = .1;   % original .01

options_.EST.UB.eta = 0.9;  
options_.EST.UB.gamma = 0.91; 
options_.EST.UB.phi = 0.96; 
options_.EST.UB.psi_N = 100; 
options_.EST.UB.rhozeta = 0.99; 
options_.EST.UB.sigmazeta = 10;  
options_.EST.UB.rho_lambda = .99;
options_.EST.UB.lambda_ss = 0.5;


%% restrictions for grid search
% options_.EST.LB.sigmazeta = 0.75;
% options_.EST.LB.L_ss = 1.5;
% options_.EST.UB.L_ss = 2.5;
% options_.EST.UB.psi_N = 50;

% that way we do one good search (when it's reasonably large)

%% Set distance (not always needed, but useful for grid search)
% No need to edit this 

options_.EST.DIST = [];
FF = options_.EST.variables;
for ii = 1:length(FF)
    var_ii = char(FF(ii));
    options_.EST.DIST.(var_ii) = options_.EST.UB.(var_ii) - options_.EST.LB.(var_ii);
end


