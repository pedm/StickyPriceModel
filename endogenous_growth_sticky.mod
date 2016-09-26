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
eta    = 0.90;       % 0.375;                   % Curvature of innovations production in R&D expenditure (original = 0.33)
gamma  = 0.0396;       % 0.35;                    % weight of current consumption on JR term; indexes strength of wealth effects (0->no wealth effect (GHH), 1-> KPR prefs)
phi    = 0.9862;       % 0.875;                   % Survival rate of technologies
% lambda = 0.4466;       % 0.075;                   % Adoption probability

% SHOCKS
rhozeta    = 0.5487; % 0.5; 
rhozeta2   = 0.0004; % 0.1;                 % Note: there's a minus in front of this (also, in estimation, must be greater than 0)
sigmazeta  = 1.0613; % 3.5;
zetabar    = 0.4777;

% NEW VARIABLES
% TODO: any suggested calibration?
gamma_pi = 1.5;
gamma_y  = 0.1;
theta    = .779;

omega    = 4.167;  
mkup_ss  = omega / (omega - 1);         % Markup. In the flex price model, markup is exogenous and given by M = ω/(ω − 1). I took this numbers from Gertler-Karadi “a model of unconventional monetary policy�?, who take them from estimates by Primiceri et al
psi_N   = 21.8893;                      % Adjustment cost to N
psi_I   = 1;                            % Adjustment cost to I

% SEPTEMBER ADDITIONS
lambda_bar = 0.5;
rho_lambda = 0.95;       % 0 < rho_lambda < 1

% Note: gg is not set here, as it depends on the steady state. 
% The param gg is instead defined in endogenous_growth_steadystate.m

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

steady;
check;

% Set seed for simulation
set_dynare_seed(092677);
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
else
    
%=========================================================================%
%%%%                       OPTIMIZATION                                %%%%
%=========================================================================%
% Much of this code comes from Bonn and Pfeifer 2014 replication files

% Starting point (based on earlier calibration)
x_start=[eta, gamma, phi, lambda, psi_N, rhozeta, rhozeta2, sigmazeta, zetabar];
x_start_unbounded = boundsINV(x_start);

% Optimizer options
H0 = 1e-2*eye(length(x_start)); % Initial Hessian 
H0 = 1e-1*eye(length(x_start)); % Initial Hessian 

crit = 1e-7; % Tolerance
nit = 1000; % Number of iterations
nit = 500;

% Make sure Dynare does not print out stuff during runs
options_.nocorr=1;
options_.noprint=1;
options_.verbosity=0;

% options_.qz_criterium = 1+1e-6; % required because it is empty by default, leading to a crash in k_order_pert
[fhat, params_unbounded] = csminwel(@distance_fcn     ,x_start_unbounded,H0,[],crit,nit);
%                          csminwel(fcn               ,x0,H0,grad,crit,nit,varargin)

%=========================================================================%
%%%%                       PLOT SOLUTION                               %%%%
%=========================================================================%

[ params ] = bounds( params_unbounded );
set_param_value('eta', params(1) );
set_param_value('gamma', params(2) );
set_param_value('phi', params(3) );
set_param_value('lambda', params(4) );
set_param_value('psi_N', params(5) );
set_param_value('rhozeta', params(6) );
set_param_value('rhozeta2', params(7) );
set_param_value('sigmazeta', params(8) );
set_param_value('zetabar', params(9) );

var_list_=[];
info = stoch_simul(var_list_);
post_processing_irfs;                                                       % Create IRFs with trend
post_processing_irfs_plot;                                                  % Plot IRFs
axis tight;
% post_processing_irfs_distance;                                              % Compute distance between model and VAR IRFs
legend('initial','var', 'b', 'b', 'final');

disp('[eta, gamma, phi, lambda, psi_N, rhozeta, rhozeta2, sigmazeta, zetabar]')
params
x_start

end