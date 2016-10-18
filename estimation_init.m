%% Initialize Estimation
options_.EST = [];

%% 1. IRF Length (you can choose between 11 and 31)
options_.EST.irf_length = 11; 

%% 2. Which parameters to estimate
options_.EST.variables = {'eta', 'phi', 'lambda_ss', 'psi_N', 'rhozeta', 'sigmazeta', 'gamma'}; 

%% 3. Define parameter bounds

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
options_.EST.UB.phi = 0.999; % original .95
options_.EST.UB.psi_N = 200; 
options_.EST.UB.rhozeta = 0.99; 
options_.EST.UB.sigmazeta = 10;  
options_.EST.UB.rho_lambda = .99;
options_.EST.UB.lambda_ss = 0.5;

%% 4. Distance between bounds
% No need to edit this

% (Only used in grid search)

options_.EST.DIST = [];
FF = options_.EST.variables;
for ii = 1:length(FF)
    var_ii = char(FF(ii));
    options_.EST.DIST.(var_ii) = options_.EST.UB.(var_ii) - options_.EST.LB.(var_ii);
end


