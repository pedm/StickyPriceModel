% Sticky Price Model of Endogenous Growth  
% August/Sept 2016 Update
% Includes Jaimovich-Rebelo Preferences, a simplified innovation sector, and adjustment costs
% Requires endogenous_growth_sticky_steadystate.m

% If you comment out this command, this code will plot the IRFs on top of your previous plots
% close all;

close all;
%===================================================================%
%                    DECLARATION OF VARIABLES                       %
%===================================================================%

do_estimate = 1; % if 0, just simulate

var

% ENDOGENOUS VARIABLES
g            ${g}$   (long_name='Growth')  % notes ...
ZD           ${Z_D}$                       % Total Measure of Technologies (adopted or not)
lambda       ${\lambda}$                     % Adoption probability
VD           ${V_D}$                       % New innovation / technologies created this period (adopted or not)
M            ${M}$                         %
J            ${J}$                         % Value of unadopted tech
H            ${H}$                         % Value of adopted tech
Pi           ${\Pi}$                       % 
ND           ${N_D}$                       % 
YDW          ${Y^W_D}$                     % New in model 4.4.1. Why is this needed? (might be needed later in the sticky price model)
YD           ${Y_D}$                       % 
CD           ${C_D}$                       % 
Lambda       ${\Lambda}$                   % Note Lambda(+1) = \Lambda_{t,t+1}
UCD          ${U_{CD}}$                    % 
muD          ${{\mu}_{D}}$                 % 
L            ${L}$                         % 
GammaD       ${\Gamma_D}$                  % 
KD           ${K_D}$                       % 
Q            ${Q}$                         % 
ID           ${I_D}$                       % 
zeta         ${\zeta}$                     % 
SD           ${\mathcal{S}_{D}}$           % 
XD           ${X_D}$                       % 
RD           ${\mathcal{R}_{D}}$           % 

% Sticky Price Variables
mkup         $\mathcal{M}$                 % Markup
pi           ${\pi}$                       % Inflation
pi_star      ${\pi^*}$                     % Inflation in the optimal reset price (or is it slighly different? in eqn 213, P_t-1 is not P*_t-1)
x1D          ${x_{1D}}$                    % Simplification 1 used in Optimal Reset Price
x2D          ${x_{2D}}$                    % Simplification 2 used in Optimal Reset Price
R_nom        ${R}$                         % R_t is nominal interest rate between t and t+1

% "Functions"
% Note: these are technically detrended
f_fcn             ${\left.       f\left( \cdot \right)            \right|}$
f_fcn_prime       ${\left.       f^‎{\prime}\left( \cdot \right)   \right|}$
g_fcn             ${\left.       g\left( \cdot \right)            \right|}$
g_fcn_prime       ${\left.       g^‎{\prime}\left( \cdot \right)   \right|}$;

% EXOGENOUS SHOCKS
varexo epsilon_chi    ${\epsilon}^{\chi}$;  

%=========================================================================%
%%%%                    DECLARATION OF PARAMETERS                      %%%%
%=========================================================================%
 
parameters 

% ORIGINAL PARAMETERS 
beta           $\beta$         
alpha          $\alpha$        
epsilon        $\epsilon$             % Inverse Frisch labor supply elasticity
rho            $\rho$                 % Inverse of intertemporal elasticity of substitution
delta          $\delta$               
chi            $\chi$                 % disutility of labor supply
vartheta       $\vartheta$              
gamma          $\gamma$               % weight of current consumption on JR term; indexes strength of wealth effects (0->no wealth effect (GHH), 1-> KPR prefs)
phi            $\phi$                 
eta            $\eta$                 
rhozeta        ${\rho}_{\zeta}$       % persistence of exogenous "innovation productivity" shock
rhozeta2       ${\rho}_{\zeta2}$      % AR(2) parameter (after a minus sign)
sigmazeta      ${\sigma}_{\zeta}$     % size of impulse on exogenous "innovation productivity" shock
zetabar        $\overline{\zeta}$

% NEW PARAMETERS (FLEX PRICE)
mkup_ss        ${\mathcal{M}^{ss}}$   % Markup
psi_N          $\psi_N$               % Magnitude of N adjustment cost
psi_I          $\psi_I$               % Magniturde if I investment cost
gg             $g^{BGP}$              % Steady state growth rate (seen in adjustment cost functions)

% NEW PARAMETERS (STICKY PRICE)
R_nom_ss       ${R^{ss}}$             % Nominal interest rate in steady state
gamma_pi       ${{\gamma}_{\pi}}$
gamma_y        ${{\gamma}_{y}}$
omega          ${\omega}$
theta          ${\theta}$

% SEPTEMBER ADDITIONS
lambda_bar     ${\bar{\lambda}}$
rho_lambda     ${\rho_\lambda}$       % 0 < rho_lambda < 1
;

%=========================================================================%
%%%%                     PARAMETERS' VALUES                            %%%%
%=========================================================================%

beta    = 0.96; 
alpha   = 0.33;
epsilon = 1/2;                   % Inverse Frisch labor supply elasticity
rho     =   1;                   % Inverse of intertemporal elasticity of substitution
% varphi  = 0.5;                 % Elasticity of Q to I/K ratio 
% varphi  = 0.5 / 10000;         % Paramater value in gensys - (albert) I think we can eliminate varphi, right?
delta   = 0.10;
chi     = 1.5652;                % Disutility of labor supply
vartheta = 1 + 1/(1-alpha);

% GROWTH PARAMETERS
eta    = 0.3;       % 0.375;                   % Curvature of innovations production in R&D expenditure (original = 0.33)
gamma  = 0.5;       % 0.35;                    % weight of current consumption on JR term; indexes strength of wealth effects (0->no wealth effect (GHH), 1-> KPR prefs)
phi    = 0.95;       % 0.875;                   % Survival rate of technologies
% lambda = 0.4466;       % 0.075;                   % Adoption probability

% SHOCKS
rhozeta    = 0.4; % 0.5; 
rhozeta2   = 0.0004; % 0.1;                 % Note: there's a minus in front of this (also, in estimation, must be greater than 0)
sigmazeta  = 0.5; % 3.5;
zetabar    = 0.6;

% NEW VARIABLES
% TODO: any suggested calibration?
gamma_pi = 1.5;
gamma_y  = 0.1;
theta    = 1;

omega    = 4.167;  
mkup_ss  = omega / (omega - 1);         % Markup. In the flex price model, markup is exogenous and given by M = ω/(ω − 1). I took this numbers from Gertler-Karadi “a model of unconventional monetary policy�?, who take them from estimates by Primiceri et al
psi_N   = 10;                      % Adjustment cost to N
psi_I   = 1;                            % Adjustment cost to I

% SEPTEMBER ADDITIONS
lambda_bar = 1.5;
rho_lambda = 0.7;       % 0 < rho_lambda < 1

% Note: gg is not set here, as it depends on the steady state. 
% The param gg is instead defined in endogenous_growth_steadystate.m

%% NEW GUESS
% eta =       0.2000;
% gamma =     0.2000;
% phi =     0.5980;
% lambda_bar =     2.0040;
% psi_N =    80.0000;
% rhozeta =     1.2000;
% rhozeta2 =     0.7920;
% sigmazeta =     2.4000;
% zetabar =     2.0800;

% eta        =    0.9530;
% gamma        =    0.7589;
% phi        =    0.9416;
% lambda_bar        =    1.5305;
% psi_N        =    15;
% rhozeta        =    0.5;
% rhozeta2        =    0.0004;
% sigmazeta        =    0.8;
% zetabar        =    0.2617;
% rho_lambda        =    0.7000;

eta            =     0.6331;
gamma            =   0.9769;
phi            =   0.9890;
lambda_bar            =   0.3;
psi_N            =  20;
rhozeta            =   0.6557;
rhozeta2            =      0.0004;
sigmazeta            =   0.9087;
rho_lambda            =   0.7905;


% Grid Search Results (1):
eta = 0.09999;
gamma = 0.899911;
phi = 0.941;
lambda_bar = 1.0045;
psi_N = 10;
rhozeta = 0.15009;
% rhozeta2 = 0.89101;
sigmazeta = 1.45;
rho_lambda = 0.892;

% Estimation results using Grid Search Results (1) as starting point:
% These params don't have a steady state... unless I modify start point in
% ss solver
% eta = 0.2602055353;
% gamma = 0.949986429;
% phi = 0.9870225118;
% lambda_bar = 0.09206466431;
% psi_N = 9.999704175;
% rhozeta = 0.8243755448;
% rhozeta2 = 0.0003999698952;
% sigmazeta = 0.4447255525;
% rho_lambda = 0.7168501364;


% Estimation Results when I include 'steady;' command in distance_fcn
% eta = 0.06983957799;
% gamma = 0.8833204851;
% phi = 0.986685405;
% lambda_bar = 0.3615669612;
% psi_N = 10.00021269;
% rhozeta = 0.6684293896;
% rhozeta2 = 0.0004;
% sigmazeta = 0.7063087393;
% rho_lambda = 0.8823283013;


% Grid Search Results:
% eta = 0.05;
% gamma = 0.9;
% phi = 0.9;
% lambda_bar = 0.2;
% % psi_N = 0;
% % rhozeta = 0;
% % rhozeta2 = 0.10005;
% % sigmazeta = 5.25;
% rho_lambda = 0.7;

% Estimation Results (Oct 5)
% eta = 0.04372746517;
% gamma = 0.9998401381;
% phi = 0.9899956566;
% lambda_bar = 0.2291965598;
% psi_N = 6.651634179;
% rhozeta = 0.8779613822;
% rhozeta2 = 0.0003455172107;
% sigmazeta = 0.3618846256;
% rho_lambda = 0.8497168674;

% TODO: why did the param results from csminwel() not have a ss?
% Shouldnt csminwel pick up on this?
% Also for some reason the plot worked...!

% Perhaps if I run GSR(1), plot it, then go to the most recent params,
% perhaps it will work then. If so, why is that? Very weird....

% TODO: try estimation again with GSR(1) as starting point. Perhaps print
% out param results with more accuracy. Does that help???


%=========================================================================%
%%%%                     EQUILIBRIUM CONDITIONS                        %%%%
%=========================================================================%
model;

%=========================ENDOGENOUS VARIABLES============================%

% 1. Evolution of number of firm varieties
% g = phi* (1 + lambda * (ZD + VD - 1));
g = lambda * phi * (ZD-1) + phi;

% 2. Evolution of technological frontier Z
% ZD * g(-1) = phi * ( ZD(-1) + (1-lambda) * VD(-1) );
ZD * g(-1) = phi * ZD(-1) + VD(-1);

% 3. Innovators' production function. Q: Why is this not the supply curve of new firms?
% VD = zetabar * zeta * ZD ^ (1-eta) * ND ^ eta ;
VD = zetabar * zeta * ZD * (ND * g(-1) / KD(-1))^eta;

% 4. Euler equation for entrepreneurs
% J =  lambda * H + (1 - lambda) * phi * Lambda(+1) * J(+1);
J = -M + phi * Lambda(+1) * (lambda*H(+1) + (1-lambda)*J(+1)   );

% 5. Price of adopted technology H
H = Pi + phi * ( Lambda(+1) * H(+1) );

% 6. Profits per period
Pi = (1/vartheta) * (1/mkup) * YDW;

% 7. Innovators' FOC wrt N
% zetabar * eta * J * zeta * ( ZD / ND )^(1-eta) =  1 + log(f_fcn_prime) *  (ND * g(-1) /  ND(-1)) +  log(f_fcn) - Lambda(+1) * log(f_fcn_prime(+1)) * (ND(+1) * g / ND )^2;
Lambda(+1) * J(+1) * zetabar * zeta * ZD * (1 /  (KD(-1)/g(-1))^eta)   * (1 / ND^(1-eta) ) =  1 + log(f_fcn_prime) *  (ND * g(-1) /  ND(-1)) +  log(f_fcn) - Lambda(+1) * log(f_fcn_prime(+1)) * (ND(+1) * g / ND )^2;

% 8.
rho_lambda * lambda_bar * phi * Lambda(+1) * (H(+1) - J(+1)) = M ^ (1 - rho_lambda);

% 9
lambda = lambda_bar * M ^ rho_lambda;

% 10. Aggregate production function
YD = YDW;

% 11. Aggregate production function W
YDW = ((KD(-1) / g(-1))^alpha) * L^(1-alpha);

% 12. Resource constraint, with adjustment cost
% YD = CD + (1 + log(g_fcn)) * ID + ND;
YD = CD + (1 + log(g_fcn)) * ID + (1 + log(f_fcn)) * ND + (ZD-1)*M;

% 13. HH's stochastic discount factor
Lambda = ((beta * UCD) / UCD(-1)) * g(-1)^(-rho);

% 14. Marginal utility of consumption
% Original:
% UCD = ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) + muD * gamma * (GammaD(-1) / ( CD * g(-1) )) ^ (1-gamma);
% To make muD positive, multiply by -1
UCD = ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) + -1*muD * gamma * (GammaD(-1) / ( CD * g(-1) )) ^ (1-gamma);

% 15. Lagrange multiplier on labor disutility law of motion (new)
% Original:
% muD   = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) - ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
% To make muD positive, replace - with +
muD   = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) + ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);

% 16. Labor market equilibrium
chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) = (1/mkup) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L);

% 17. Labor disutility term
GammaD = (CD^gamma) * (GammaD(-1) / g(-1) )^(1-gamma);

% 18. Capital evolution
KD = (1-delta) * (KD(-1) / g(-1)) + ID ;

% 19. Q-equation (capital producers)
Q = 1 + log(g_fcn) + ((ID * g(-1)) / ID(-1)) * log(g_fcn_prime) - Lambda(+1) * ((ID(+1) * g) / ID)^2 * log(g_fcn_prime(+1));

% 20. Equation 263. Capital Euler Equation
Q = Lambda(+1) * ((g* (vartheta - 1) *YDW(+1) * alpha)/(mkup * KD * vartheta) + Q(+1) * (1 - delta));

% 21. Exogenous shock to entrepreneurs' production function
log(zeta) = rhozeta * log(zeta(-1)) - rhozeta2 * log(zeta(-2)) + sigmazeta * epsilon_chi;

% 22. Stock market value (taken directly from May model)
SD = Q * KD + H  +  J * ( ZD + VD - 1 )  + XD;

% 23. 
XD =  ( Lambda(+1) * g * ( J(+1) * VD(+1) + XD(+1) ) );

% 24. Aggregate R&D expenditure [ (albert) in this model, R&D expenditures
% are given by the variable N)
RD = ND;   

%====================== ADJUSTMENT COST FCNS =============================%

f_fcn = exp( (psi_N / 2) * ((ND * g(-1) / ND(-1) ) - gg)^2 );
f_fcn_prime = exp( psi_N * ((ND * g(-1) / ND(-1) ) - gg)   );

g_fcn = exp( (psi_I / 2) * ((ID * g(-1) / ID(-1) ) - gg)^2 );
g_fcn_prime = exp( psi_I * ((ID * g(-1) / ID(-1) ) - gg)   );

%=========================STICKY PRICE EQNS===============================%

% Pricing Equation
pi ^(1-omega) = theta + (1-theta)*pi_star^(1-omega);

% Optimal Pricing Equation 
pi_star = (omega / (omega - 1)) * (x1D / x2D) * pi;

% Simplifications used in the Optimal Pricing Equation
x1D = UCD * (1/mkup) * YD + beta * theta * g^(1-rho) * pi(+1)^omega * x1D(+1);

x2D = UCD * YD + beta * theta * g^(1-rho) * pi(+1)^(omega-1) * x2D(+1);

% Euler Equation
1 = Lambda(+1) * R_nom / pi(+1);

% Taylor Rule
R_nom / R_nom_ss = pi ^ gamma_pi * (  (1/mkup) / (1/mkup_ss)   )^gamma_y;

end; 

write_latex_dynamic_model;
write_latex_static_model;

%=========================================================================%
%%%%                       RUN                                         %%%%
%=========================================================================%

% Load PVAR IRFs
global pvarcoirfs_clean;
load 'pvar_coirfs_full';
pvarcoirfs_clean = pvarcoirfs;

shocks;
var epsilon_chi;
stderr 1;
end;

% Setup counters
COUNT.total = 0;
COUNT.ss_found = 0;
COUNT.ss_notfound = 0;
COUNT.ss_neg = 0;

% Set seed for simulation
set_dynare_seed(092677);

evalin('base','save level0workspace oo_ M_ options_')

steady;
check;

% Produce simulation using above calibration, compare with VAR IRFs
% NOTE: loglinear option causes oo_.steady_state to become logged
stoch_simul(order=1,periods=600, irf=11, nograph, nodisplay, nocorr, nomoments, loglinear);
post_processing_irfs;                                                       % Create IRFs with trend
plot_var_irfs;                                                              % Plot VAR IRFs
post_processing_irfs_plot;                                                  % Plot IRFs
post_processing_irfs_distance;                                              % Compute distance between model and VAR IRFs


% Change parameters, solve again, and plot
% set_param_value('zetabar', 0.6);
% stoch_simul(order=1,periods=600, irf=11, nograph, nodisplay, nocorr, nofunctions, nomoments, noprint, loglinear);
% post_processing_irfs;                                                       % Create IRFs with trend
% post_processing_irfs_plot;                                                  % Plot IRFs
% post_processing_irfs_distance;                                              % Compute distance between model and VAR IRFs
% 
% % Change parameters, solve again, and plot
% set_param_value('sigmazeta', 5);
% stoch_simul(order=1,periods=600, irf=11, nograph, nodisplay, nocorr, nofunctions, nomoments, noprint, loglinear);
% post_processing_irfs;                                                       % Create IRFs with trend
% post_processing_irfs_plot;                                                  % Plot IRFs
% post_processing_irfs_distance;                                              % Compute distance between model and VAR IRFs



% You can copy and paste the above lines in order to continue playing around with the calibration
% set_param_value('eta', 0.3);
% stoch_simul(order=1,periods=600, irf=11, nograph, nodisplay, nocorr, nofunctions, nomoments, noprint, loglinear);
% post_processing_irfs;                                                       % Create IRFs with trend
% post_processing_irfs_plot;                                                  % Plot IRFs
% post_processing_irfs_distance;                                              % Compute distance between model and VAR IRFs

if do_estimate == 0;
    return;
elseif do_estimate == 1
    
%=========================================================================%
%%%%                       OPTIMIZATION                                %%%%
%=========================================================================%
% Much of this code comes from Bonn and Pfeifer 2014 replication files

% Starting point (based on earlier calibration)
x_start=[eta, gamma, phi, lambda_bar, psi_N, rhozeta, rhozeta2, sigmazeta, rho_lambda]; % zetabar, 
x_start_unbounded = boundsINV(x_start);

% Optimizer options
H0 = 1e-2*eye(length(x_start)); % Initial Hessian 
H0 = 1e-1*eye(length(x_start)); % Initial Hessian 

crit = 1e-7; % Tolerance
nit = 1000; % Number of iterations
nit = 800;

% nit = 20;
nit = 200;

% Make sure Dynare does not print out stuff during runs
options_.nocorr=1;
options_.noprint=1;
options_.verbosity=0;

% options_.qz_criterium = 1+1e-6; % required because it is empty by default, leading to a crash in k_order_pert
[fhat, params_unbounded] = csminwel(@distance_fcn     ,x_start_unbounded,H0,[],crit,nit);
%                          csminwel(fcn               ,x0,H0,grad,crit,nit,varargin)

fhat

%=========================================================================%
%%%%                       ALT ESTIMATION                              %%%%
%=========================================================================%
elseif do_estimate == 2
    % Starting point (based on earlier calibration)
    x_start=[eta, gamma, phi, lambda_bar, psi_N, rhozeta, rhozeta2, sigmazeta, rho_lambda]; % zetabar, 
    x_start_unbounded = boundsINV(x_start);
    
    opts = psoptimset('Display','diagnose'); % debugging % , 'MaxIter', 20
    % opts = psoptimset('Display', 'off'); % estimation
    [params_unbounded, FVAL,EXITFLAG] = patternsearch(@distance_fcn, x_start_unbounded,[],[],[],[],[],[],[],opts)
end

%=========================================================================%
%%%%                       PLOT SOLUTION                               %%%%
%=========================================================================%
if do_estimate > 0
[ params ] = bounds( params_unbounded );
set_param_value('eta', params(1) );
set_param_value('gamma', params(2) );
set_param_value('phi', params(3) );
set_param_value('lambda_bar', params(4) );
set_param_value('psi_N', params(5) );
set_param_value('rhozeta', params(6) );
set_param_value('rhozeta2', params(7) );
set_param_value('sigmazeta', params(8) );
% set_param_value('zetabar', params(9) );
set_param_value('rho_lambda', params(9) );


oo_.irfs = {}; % erase the old IRFs
steady;
var_list_=[];
info = stoch_simul(var_list_);                                              % WARNING: this does not compute the steady state. It just uses the pre defined ss
post_processing_irfs;                                                       % Create IRFs with trend
post_processing_irfs_plot;                                                  % Plot IRFs
axis tight;
% post_processing_irfs_distance;                                              % Compute distance between model and VAR IRFs
legend('initial','var', 'b', 'b', 'final');

disp('[eta, gamma, phi, lambda_bar, psi_N, rhozeta, rhozeta2, sigmazeta, zetabar]')
params
x_start

% Print Estimation Results
disp(sprintf('eta = %0.10g;', params(1) ));
disp(sprintf('gamma = %0.10g;', params(2) ));
disp(sprintf('phi = %0.10g;', params(3) ));
disp(sprintf('lambda_bar = %0.10g;', params(4) ));
disp(sprintf('psi_N = %0.10g;', params(5) ));
disp(sprintf('rhozeta = %0.10g;', params(6) ));
disp(sprintf('rhozeta2 = %0.10g;', params(7) ));
disp(sprintf('sigmazeta = %0.10g;', params(8) ));
disp(sprintf('rho_lambda = %0.10g;', params(9) ));
end
