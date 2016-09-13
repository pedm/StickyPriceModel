% Flex Price Model of Endogenous Growth
% August/Sept 2016 Update
% Includes Jaimovich-Rebelo Preferences, a simplified innovation sector, and adjustment costs
% Requires endogenous_growth_steadystate.m

% If you comment out this command, this code will plot the IRFs on top of your previous plots
% close all;

close all;
%===================================================================%
%                    DECLARATION OF VARIABLES                       %
%===================================================================%

do_estimate = 1; % if 0, just simulate

var
%% comment
% ENDOGENOUS VARIABLES
g            ${g}$                         % notes ...
ZD           ${Z_D}$                       % Total Measure of Technologies (adopted or not)
VD           ${V_D}$                       % New innovation / technologies created this period (adopted or not)
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
% varphi         $\varphi$              % elasticity of Q to I/K ratio 
delta          $\delta$               
chi            $\chi$                 % disutility of labor supply
vartheta       $\vartheta$              
gamma          $\gamma$               % weight of current consumption on JR term; indexes strength of wealth effects (0->no wealth effect (GHH), 1-> KPR prefs)
phi            $\phi$                 
eta            $\eta$                 
lambda         $\lambda$              % adoption probability
rhozeta        ${\rho}_{\zeta}$       % persistence of exogenous "innovation productivity" shock
rhozeta2       ${\rho}_{\zeta2}$      % AR(2) parameter (after a minus sign)
sigmazeta      ${\sigma}_{\zeta}$     % size of impulse on exogenous "innovation productivity" shock
zetabar        $\overline{\zeta}$

% NEW PARAMETERS
M              $\mathcal{M}$          % Markup
psi_N          $\psi_N$               % Magnitude of N adjustment cost
psi_I          $\psi_I$               % Magniturde if I investment cost
gg             $g^{BGP}$              % Steady state growth rate (seen in adjustment cost functions)
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
eta    = 0.999;       % 0.375;                   % Curvature of innovations production in R&D expenditure (original = 0.33)
gamma  = 0.0396;       % 0.35;                    % weight of current consumption on JR term; indexes strength of wealth effects (0->no wealth effect (GHH), 1-> KPR prefs)
phi    = 0.9862;       % 0.875;                   % Survival rate of technologies
lambda = 0.4466;       % 0.075;                   % Adoption probability

% SHOCKS
rhozeta    = 0.5487; % 0.5; 
rhozeta2   = 0.0004; % 0.1;                 % Note: there's a minus in front of this (also, in estimation, must be greater than 0)
sigmazeta  = 1.0613; % 3.5;
zetabar    = 0.4777;


% NEW VARIABLES
M       = 4.167 / (4.167 - 1);         % Markup. In the flex price model, markup is exogenous and given by M = ω/(ω − 1). I took this numbers from Gertler-Karadi “a model of unconventional monetary policy�?, who take them from estimates by Primiceri et al
psi_N   = 21.8893;                           % Adjustment cost to N
psi_I   = 1;                           % Adjustment cost to I

% Note: gg is not set here, as it depends on the steady state. 
% The param gg is instead defined in endogenous_growth_steadystate.m

%=========================================================================%
%%%%                     EQUILIBRIUM CONDITIONS                        %%%%
%=========================================================================%
model;

%=========================ENDOGENOUS VARIABLES============================%

% 1. Evolution of number of firm varieties
g = phi* (1 + lambda * (ZD + VD - 1));

% 2. Evolution of technological fronteir Z
ZD * g(-1) = phi * ( ZD(-1) + (1-lambda) * VD(-1) );

% 3. Innovators' production function. Q: Why is this not the supply curve of new firms?
% PAT'S MODIFICATION (I replace exp(zeta) with zeta to account for zeta being 1 in steady state)
VD = zetabar * zeta * ZD ^ (1-eta) * ND ^ eta ;

% 4. Euler equation for entrepreneurs
J =  lambda * H + (1 - lambda) * phi * Lambda(+1) * J(+1);

% 5. Price of adopted technology H
H = Pi + phi * ( Lambda(+1) * H(+1) );

% 6. Profits per period
% NOTE: in the flex price model, M is an endogenous var
Pi = (1/vartheta) * (1/M) * YDW;

% 7. Innovators' FOC wrt N
% PAT'S MODIFICATION (I replace exp(zeta) with zeta to account for zeta being 1 in steady state)
% log() added so that f_fcn can be 1 in ss
zetabar * eta * J * zeta * ( ZD / ND )^(1-eta) =  1 + log(f_fcn_prime) *  (ND * g(-1) /  ND(-1)) +  log(f_fcn) - Lambda(+1) * log(f_fcn_prime(+1)) * (ND(+1) * g / ND )^2;

% 8. Aggregate production function
YD = YDW;

% 9. Aggregate production function W
YDW = ((KD(-1) / g(-1))^alpha) * L^(1-alpha);

% 10. Resource constraint, with adjustment cost
% log() added so that g_fcn can be 1 in ss
YD = CD + (1 + log(g_fcn)) * ID + ND;

% 11. HH's stochastic discount factor
Lambda = ((beta * UCD) / UCD(-1)) * g(-1)^(-rho);

% 12. Marginal utility of consumption
% Make muD positive, multiply by -1
UCD = ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) + -1*muD * gamma * (GammaD(-1) / ( CD * g(-1) )) ^ (1-gamma);

% 13. Lagrange multiplier on labor disutility law of motion (new) (equation 258)
% Albert: equation 258 needs to be corrected - it will need a “g�? term as well.
% To use dynare notation, here's the change C(+1) / Gamma = CD(+1) * g / GammaD
% Original:
% muD = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) - ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
% New so that muD is positive
muD   = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) + ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);

% 14. Labor market equilibrium (eqn 259)
% NOTE: in this equation, it seems that M is an exogenous param
chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) = (1/M) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L);

% 15. Labor disutility term
GammaD = (CD^gamma) * (GammaD(-1) / g(-1) )^(1-gamma);

% 16. Capital evolution
KD = (1-delta) * (KD(-1) / g(-1)) + ID ;

% 17. Q-equation (capital producers)
% log() added so that g_fcn can be 1 in ss
Q = 1 + log(g_fcn) + ((ID * g(-1)) / ID(-1)) * log(g_fcn_prime) - Lambda(+1) * ((ID(+1) * g) / ID)^2 * log(g_fcn_prime(+1));

% PREVIOUS: 17. Q-equation
% Q = ((ID * g(-1) * (I_K_ss^-1)) / (KD(-1) )  )^varphi;

% 18. Equation 263. Capital Euler Equation (perhaps?)
% NOTE: in this equation, it seems that M is an exogenous param
Q = Lambda(+1) * ((g* (vartheta - 1) *YDW(+1) * alpha)/(M * KD * vartheta) + Q(+1) * (1 - delta));

% 19. Exogenous shock to entrepreneurs' production function
log(zeta) = rhozeta * log(zeta(-1)) - rhozeta2 * log(zeta(-2)) + sigmazeta * epsilon_chi;

% 20. Stock market value (taken directly from May model)
% original
% SD =  Q * KD + (H - Pi)  +  J * ( ZD - 1 )  + XD;
% new
SD = Q * KD + H  +  J * ( ZD + VD - 1 )  + XD;

% 21. 
XD =  ( Lambda(+1) * g * ( J(+1) * VD(+1) + XD(+1) ) );

% 22. Aggregate R&D expenditure [ (albert) in this model, R&D expenditures
% are given by the variable N)
RD = ND;   

% ADJUSTMENT COST FUNCTIONS
f_fcn = exp( (psi_N / 2) * ((ND * g(-1) / ND(-1) ) - gg)^2 );
f_fcn_prime = exp( psi_N * ((ND * g(-1) / ND(-1) ) - gg)   );

g_fcn = exp( (psi_I / 2) * ((ID * g(-1) / ID(-1) ) - gg)^2 );
g_fcn_prime = exp( psi_I * ((ID * g(-1) / ID(-1) ) - gg)   );

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
return
check;

% Set seed for simulation
set_dynare_seed(092677);
% Produce simulation using above calibration, compare with VAR IRFs
% NOTE: loglinear option causes oo_.steady_state to become logged
stoch_simul(order=1,periods=600, irf=11, nograph, nodisplay, nocorr, nomoments, loglinear);
post_processing_irfs;                                                       % Create IRFs with trend
post_processing_irfs_plot;                                                  % Plot IRFs
post_processing_irfs_distance;                                              % Compute distance between model and VAR IRFs
plot_var_irfs;                                                              % Plot VAR IRFs


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

stoch_simul(order=1,periods=600, irf=11, nograph, nodisplay, nocorr, nofunctions, nomoments, noprint, loglinear);
post_processing_irfs;                                                       % Create IRFs with trend
post_processing_irfs_plot;                                                  % Plot IRFs
axis tight;
% post_processing_irfs_distance;                                              % Compute distance between model and VAR IRFs
legend('initial','var', 'b', 'b', 'final');

disp('[eta, gamma, phi, lambda, psi_N, rhozeta, rhozeta2, sigmazeta, zetabar]')
params
x_start

end