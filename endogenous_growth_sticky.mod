%% Sticky Price Model of Endogenous Growth  
% August/Sept 2016 Update
% Includes Jaimovich-Rebelo Preferences, a simplified innovation sector, and adjustment costs
% Requires endogenous_growth_sticky_steadystate.m
% Albert SS Solver Oct 2016

close all;
%===================================================================%
%                    DECLARATION OF VARIABLES                       %
%===================================================================%
var

% ENDOGENOUS VARIABLES
g            ${g}$   (long_name='Growth')  % 
ZD           ${Z_D}$                       % Total Measure of Technologies (adopted or not)
lambda       ${\lambda}$                   % Adoption probability
VD           ${V_D}$                       % New innovation / technologies created this period (adopted or not)
M            ${M}$                         %
J            ${J}$                         % Value of unadopted tech
H            ${H}$                         % Value of adopted tech
PI           ${\Pi}$                       % 
ND           ${N_D}$                       % 
YDW          ${Y^W_D}$                     % (Might be needed later in the sticky price model)
YD           ${Y_D}$                       % 
CD           ${C_D}$                       % 
LAMBDA       ${\Lambda}$                   % Household Stochastic Discount Factor (Note Lambda(+1) = \Lambda_{t,t+1} )
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
DELTAN       ${\Delta^N}$                  % financial shock affecting innovators
DELTAM       ${\Delta^M}$                  % financial shock affecting adopters

% Sticky Price Variables
mkup         ${\mathcal{M}}$               % Markup
pi           ${\pi}$                       % Inflation
pi_star      ${\pi^*}$                     % Inflation in the optimal reset price (or is it slighly different? in eqn 213, P_t-1 is not P*_t-1)
x1D          ${x_{1D}}$                    % Simplification 1 used in Optimal Reset Price
x2D          ${x_{2D}}$                    % Simplification 2 used in Optimal Reset Price
R_nom        ${R}$                         % R_t is nominal interest rate between t and t+1

% "Functions"
% Note: these are technically detrended
f_fcn             ${\left.       f\left( \cdot \right)            \right|}$
f_fcn_prime       ${\left.       f^{\prime}\left( \cdot \right)   \right|}$
g_fcn             ${\left.       g\left( \cdot \right)            \right|}$
g_fcn_prime       ${\left.       g^{\prime}\left( \cdot \right)   \right|}$;

% EXOGENOUS SHOCKS
varexo
epsilon_chi  ${\epsilon}^{\chi}$
epsilon_n    ${\epsilon}^n$
;  

%=========================================================================%
%%%%                    DECLARATION OF PARAMETERS                      %%%%
%=========================================================================%
 
parameters 

% ORIGINAL PARAMETERS 
beta           $\beta$         
alpha          $\alpha$        
epsilon        $\epsilon$             		% Inverse Frisch labor supply elasticity
rho            $\rho$                 		% Inverse of intertemporal elasticity of substitution
delta          $\delta$               		%
chi            $\chi$                 		% disutility of labor supply
vartheta       $\vartheta$              	% 
gamma          $\gamma$               		% weight of current consumption on JR term; indexes strength of wealth effects (0->no wealth effect (GHH), 1-> KPR prefs)
phi            $\phi$                 		% 
eta            $\eta$                 		% 
rhozeta        ${{\rho}_{\zeta}}$       	% persistence of exogenous "innovation productivity" shock
sigmazeta      ${{\sigma}_{\zeta}}$     	% size of impulse on exogenous "innovation productivity" shock
zeta_bar       ${\overline{\zeta}}$			% 
rhon           ${\rho}_n$             		% persistence of financial shock on innovators
sigman         ${\sigma}_n$           		% sd of financial shock on innovators
psi_N          ${\psi_N}$               	% Magnitude of N adjustment cost
psi_I          ${\psi_I}$               	% Magniturde if I investment cost
alpha_N        $\alpha_N$             		% financial spillover from innovation to adoption

% STICKY PRICE PARAMETERS
gamma_pi       ${{\gamma}_{\pi}}$
gamma_y        ${{\gamma}_{y}}$
gamma_r        ${{\gamma}_{r}}$ 
omega          ${\omega}$
theta          ${\theta}$

% SEPTEMBER ADDITIONS
lambda_bar     ${\bar{\lambda}}$			% 
rho_lambda     ${\rho_\lambda}$       		% 0 < rho_lambda < 1

% CALIBRATED STEADY-STATE VALUES
g_ss           ${g^{SS}}$              		% Steady-state growth rate (seen in adjustment cost functions)
L_ss           ${L^{SS}}$              		% Steady-state labor
lambda_ss      ${\lambda^{SS}}$        		% Steady-state adoption rate
mkup_ss        ${\mkup^{SS}}$				% Steady-state markup
;

%=========================================================================%
%%%%                     PARAMETERS' VALUES                            %%%%
%=========================================================================%

% CONVENTIONAL
beta    = 0.96;       						% discount factor
alpha   = 0.33;       						% capital share
epsilon = 0.5;        						% inv Frisch
rho     = 1;          						% IES
delta   = 0.1;        						% capital deprec

% GROWTH
vartheta  = 2.4925; 						% elast. of subst. across varieties (mkup = vartheta/(vartheta-1)
gamma     = 1.00;   		                % (Before .9) Jaimovich-Rebelo term on preferences (0-> no wealth effect)
phi       = 0.90;   						% Survival rate of technologies
eta       = 0.25;   						% elasticity of innovations to R&D
rho_lambda= 0.50;   						% elasticity of adoption prob. to adoption expenditure

% STICKY PRICE VARIABLES
gamma_pi = 1.5;
gamma_y  = 0.5;
gamma_r  =  0.32;  
theta    = 0.65;
omega    = 4.167;      		            	%  Markup. In the flex price model, markup is exogenous and given by M = ?/(? ? 1). I took this numbers from Gertler-Karadi “a model of unconventional monetary policy??, who take them from estimates by Primiceri et al

psi_I    = 1;                      			%  Adjustment cost to Inv
psi_N    = 10;                     			%  Adjustment cost to N

% CALIBRATED STEADY-STATE VALUES
g_ss      = 1.0118; 						% calibrate g (growth rate); back out zeta_bar
lambda_ss = 0.15;   						% calibrate lambda (adoption probability); back out lambda_bar
L_ss      =  1;        						% calibrate SS L; back out chi 
mkup_ss   = omega/(omega-1);

% AR(1) parameters for zeta shock (productivity of R&D)
rhozeta    = 0.90;
sigmazeta  = 0.50;

% AR(1) parameters for DeltaN shock (innovators/adopters wedge)
rhon = 0.70;
sigman = 1.00;

%% Best estimate so far
%  Will need to update, since the model and calibrated params have changed
eta = 0.08303187684;
alpha_N = 0.009623755187;
psi_N = 25.08490715;
rhon = 0.8044297245;

%=========================================================================%
%%%%                     EQUILIBRIUM CONDITIONS                        %%%%
%=========================================================================%

model;
%=========================ENDOGENOUS VARIABLES============================%
% 1. Evolution of number of firm varieties
g = lambda * phi * (ZD-1) + phi;

% 2. Evolution of technological frontier Z
ZD * g(-1) = phi * ZD(-1) + VD(-1);

% 3. Innovators' production function.
VD = zeta_bar * zeta * ZD * (ND * g(-1) / KD(-1))^eta;

% 4. Euler equation for entrepreneurs
J = -M + phi * LAMBDA(+1) * (lambda*H(+1) + (1-lambda)*J(+1)   );

% 5. Price of adopted technology H
H = PI + phi * ( LAMBDA(+1) * H(+1) );

% 6. Profits per period
PI = (1/vartheta) * (1/mkup) * YDW;

% 7. Innovators' FOC wrt N
LAMBDA(+1) * J(+1) * zeta_bar * zeta * ZD * (1 /  (KD(-1)/g(-1))^eta)   * (1 / ND^(1-eta) ) * (1 / DELTAN) =  1 + log(f_fcn_prime) *  (ND * g(-1) /  ND(-1)) +  log(f_fcn) - LAMBDA(+1) * log(f_fcn_prime(+1)) * (ND(+1) * g / ND )^2;

% 8. Adopters' FOC
rho_lambda * lambda_bar * phi * LAMBDA(+1) * (H(+1) - J(+1)) * (1 / DELTAM ) = M ^ (1 - rho_lambda);

% 9. Adoption Probability:
lambda = lambda_bar * M ^ rho_lambda;

% 10. Aggregate production function
YD = YDW;

% 11. Aggregate production function W
YDW = ((KD(-1) / g(-1))^alpha) * L^(1-alpha);

% 12. Resource constraint, with adjustment cost
YD = CD + (1 + log(g_fcn)) * ID + (1 + log(f_fcn)) * ND + (ZD-1)*M;

% 13. HH's stochastic discount factor
LAMBDA = ((beta * UCD) / UCD(-1)) * g(-1)^(-rho);

% 14. Marginal utility of consumption
% Modifed such that muD is positive
UCD = ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) - muD * gamma * (GammaD(-1) / ( CD * g(-1) )) ^ (1-gamma);

% 15. Lagrange multiplier on labor disutility law of motion (new)
% Modifed such that muD is positive
muD   = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) + ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);

% 16. Labor market equilibrium
chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) = (1/mkup) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L);

% 17. Labor disutility term
GammaD = (CD^gamma) * ( GammaD(-1) / g(-1) )^(1-gamma);

% 18. Capital evolution
KD = (1-delta) * (KD(-1) / g(-1)) + ID ;

% 19. Q-equation (capital producers)
Q = 1 + log(g_fcn) + ((ID * g(-1)) / ID(-1)) * log(g_fcn_prime) - LAMBDA(+1) * ((ID(+1) * g) / ID)^2 * log(g_fcn_prime(+1));

% 20. Equation 263. Capital Euler Equation
Q = LAMBDA(+1) * ((g* (vartheta - 1) *YDW(+1) * alpha)/(mkup * KD * vartheta) + Q(+1) * (1 - delta));

% 21. Exogenous shock to entrepreneurs' production function
log(zeta) = rhozeta * log(zeta(-1)) + sigmazeta * epsilon_chi;

% 22. Stock market value (taken directly from May model)
SD = Q * KD + H  +  J * ( ZD + VD - 1 )  + XD;

% 23. 
XD =  ( LAMBDA(+1) * g * ( J(+1) * VD(+1) + XD(+1) ) );
 
% 24. Aggregate R&D expenditure
RD = ND;   

% 25. Exogenous financial shock on innovators
log(DELTAN) = rhon * log(DELTAN(-1)) - sigman * epsilon_n;

% 26. Financial disturbance on adopters:
log(DELTAM) = alpha_N * log(DELTAN);
 
%====================== ADJUSTMENT COST FCNS =============================%
f_fcn = exp( (psi_N / 2) * ((ND * g(-1) / ND(-1) ) - g_ss)^2 ); % 22
f_fcn_prime = exp( psi_N * ((ND * g(-1) / ND(-1) ) - g_ss)   ); % 23
g_fcn = exp( (psi_I / 2) * ((ID * g(-1) / ID(-1) ) - g_ss)^2 ); % 24
g_fcn_prime = exp( psi_I * ((ID * g(-1) / ID(-1) ) - g_ss)   ); % 25

%=========================STICKY PRICE EQNS===============================%

% 26. Pricing Equation
pi ^(1-omega) = theta + (1-theta)*pi_star^(1-omega);

% 27. Optimal Pricing Equation 
pi_star = (omega / (omega - 1)) * (x1D / x2D) * pi;

% 28, 29. Simplifications used in the Optimal Pricing Equation
x1D = UCD * (1/mkup) * YD + beta * theta * g^(1-rho) * pi(+1)^omega * x1D(+1);
x2D = UCD * YD + beta * theta * g^(1-rho) * pi(+1)^(omega-1) * x2D(+1);

% 30. Euler Equation
1 = LAMBDA(+1) * R_nom / pi(+1);

% 31. Taylor Rule
R_nom  = ( R_nom(-1) )^gamma_r * ( pi ^ gamma_pi * (  (1/mkup) / (1/mkup_ss)   )^gamma_y *  1/(beta*g_ss^(-rho))  )^(1-gamma_r);

end;

% write_latex_dynamic_model;
% write_latex_static_model;

%=========================================================================%
%%%%                       RUN                                         %%%%
%=========================================================================%
shocks;
var epsilon_chi = 1;
var epsilon_n   = 1;
end;

% This saves the oo_ and M_ objects without any steady state defined (used in estimation)
evalin('base','save level0workspace oo_ M_ options_')

steady;
check;

% Produce simulation
stoch_simul(order=1,periods=600, irf=31, nograph, nodisplay, nocorr, nomoments, loglinear, irf_shocks = (epsilon_n) );

% Plot IRFs and comapre with PVAR IRFs
addpath('Scripts');
plot_var_irfs;                                                              % Plot VAR IRFs
post_processing_irfs;                                                       % Create IRFs with trend
post_processing_irfs_plot;                                                  % Plot IRFs
steady_state_targets
