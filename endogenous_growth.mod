% Flex Price Model of Endogenous Growth
% August 2016 Update
% Includes Jaimovich-Rebelo Preferences, a simplified innovation sector, and adjustment costs

%===================================================================%
%                    DECLARATION OF VARIABLES                       %
%===================================================================%

var

% ENDOGENOUS VARIABLES
g            ${g}$                         % notes ...
ZD           ${Z_D}$                       % Total Measure of Technologies (adopted or not)
VD           ${V_D}$                       % New innovation / technologies created this period (adopted or not)
J            ${J}$                         % Value of unadopted tech
H            ${H}$                         % Value of adopted tech
Pi           ${\Pi}$                       % 
ND           ${N_D}$                       % 
YDW          ${Y^W_D}$                     % New in model 4.4.1. Why is this needed?
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
RD           ${\mathcal{R}_{D}}$           % 

% "Functions"
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
varphi         $\varphi$              % elasticity of Q to I/K ratio 
delta          $\delta$               
chi            $\chi$                 % disutility of labor supply
vartheta       $\vartheta$              
gamma          $\gamma$               % elasticity of labor disutility to technology
phi            $\phi$                 
eta            $\eta$                 
LS             $LS$                   
lambda         $\lambda$              % adoption probability
rhozeta        ${\rho}_{\zeta}$       % persistence of exogenous "innovation productivity" shock

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
varphi  = 0.5;                   % elasticity of Q to I/K ratio 
varphi = 0.5 / 10000;            % param in gensys
delta   = 0.10;
chi     = 0.2236;                % disutility of labor supply
vartheta = 1 + 1/(1-alpha) ;

% GROWTH PARAMETERS
gamma  = 0.10;                   % elasticity of labor disutility to technology
gamma  = 0.12;                   % param in gensys
phi    = 0.975;
phi    = 0.90;                   % param in gensys
eta    = 0.33;
LS     = 0.0645; 
lambda = 0.1;                    % adoption probability

% SHOCKS
rhozeta = 0.5; 

% NEW VARIABLES
M       = 4.167 / (4.167 - 1);         % Flex price model: the markup M_t is exogenous and given by M = ω/(ω − 1). I took this numbers from Gertler-Karadi “a model of unconventional monetary policy”, who take them from estimates by Primiceri et al
psi_N   = 1;
psi_I   = 1;

%=========================================================================%
%%%%                   CUSTOM STEADY STATE SOLVER                      %%%%
%=========================================================================%

% Solve for steady state using these parameters
delete('ss_flexprice.mat');
disp('Steady State Solution:')
ss_solver_flex_price_model( M_.param_names, M_.params);

% Set gg param equal to g in steady state
load ss_flexprice.mat
gg = steady_g;

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
VD = zeta * ZD ^ (1-eta) * ND ^ eta;

% PREVIOUS: 3. Supply curve of new firms
% J = ((1/zeta)^(1 / eta) * (1/eta) *(1 /LS )^((1 - eta)/eta)*(VD/ZD)^((1 - eta)/eta));

% 4. Euler equation for entrepreneurs
J =  lambda * H + (1 - lambda) * phi * Lambda(+1) * J(+1);

% 5. Price of adopted technology H
H = Pi + phi * ( Lambda(+1) * H(+1) );

% 6. Profits per period
Pi = (1/vartheta) * (1/M) * YDW;

% 7. Innovators' FOC wrt N
% PAT'S MODIFICATION (I replace exp(zeta) with zeta to account for zeta being 1 in steady state)
eta * J * zeta * ( ZD / ND )^(1-eta) =  1 + f_fcn_prime *  (ND * g(-1) /  ND(-1)) +  f_fcn - Lambda(+1) * f_fcn_prime(+1) * (ND(+1) * g / ND )^2;

% PREVIOUS: 7. Final output used by the entrepreneurial sector
% ND = eta * VD * J;

% 8. Aggregate production function
YD = YDW;

% 9. Aggregate production function W
YDW = ((KD(-1) / g(-1))^alpha) * L^(1-alpha);

% 10. Resource constraint, with adjustment cost
YD = CD + (1 + g_fcn) * ID + ND;

% 11. HH's stochastic discount factor
Lambda = ((beta * UCD) / UCD(-1)) * g(-1)^(-rho);

% 12. Marginal utility of consumption
UCD = ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) + muD * gamma * (GammaD(-1) / ( CD * g(-1) )) ^ (1-gamma);

% 13. Lagrange multiplier on labor disutility law of motion (new) (equation 258)
% Albert: equation 258 needs to be corrected - it will need a “g” term as well.
% To use dynare notation, here's the change C(+1) / Gamma = CD(+1) * g / GammaD
muD = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) - ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);

% 14. Labor market equilibrium (eqn 259)
chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) = (1/M) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L);

% 15. Labor disutility term
GammaD = (CD^gamma) * (GammaD(-1) / g(-1) )^(1-gamma);

% 16. Capital evolution
KD = (1-delta) * (KD(-1) / g(-1)) + ID ;

% 17. Q-equation (capital producers)
Q = 1 + g_fcn + ((ID * g(-1)) / ID(-1)) * g_fcn_prime - Lambda(+1) * ((ID(+1) * g) / ID)^2 * g_fcn_prime(+1);

% PREVIOUS: 17. Q-equation
% Q = ((ID * g(-1) * (I_K_ss^-1)) / (KD(-1) )  )^varphi;

% 18. Equation 263
Q = Lambda(+1) * ((g* (vartheta - 1) *YDW(+1) * alpha)/(M * KD * vartheta) + Q(+1) * (1 - delta));

% 19. Exogenous shock to entrepreneurs' production function
log(zeta) = rhozeta * log(zeta(-1)) + 0.1 * epsilon_chi;

% 20. Aggregate R&D expenditure
RD = J * VD;

% ADJUSTMENT COST FUNCTIONS
f_fcn = (psi_N / 2) * ((ND * g(-1) / ND(-1) ) - gg)^2;
f_fcn_prime = psi_N * ((ND * g(-1) / ND(-1) ) - gg);

g_fcn = (psi_I / 2) * ((ID * g(-1) / ID(-1) ) - gg)^2;
g_fcn_prime = psi_I * ((ID * g(-1) / ID(-1) ) - gg);

end; 

write_latex_dynamic_model;
write_latex_static_model;

%===================================================================%
%%%%             INITIAL CONDITIONS FOR STEADY STATE             %%%%
%===================================================================%

initval;

	% Set initial guess equal to results from ss solver
	g = steady_g;
	VD  = steady_VD;
	ZD  = steady_ZD;         
	zeta  = steady_zeta;          
	ND  = steady_ND;          
	Lambda  = steady_Lambda;          
	J  = steady_J;          
	H  = steady_H;         
	Pi  = steady_Pi;          
	YDW  = steady_YDW;          
	YD  = steady_YD;          
	Q  = steady_Q;          
	KD  = steady_KD;          
	L  = steady_L;          
	ID  = steady_ID;          
	CD  = steady_CD;          
	GammaD  = steady_GammaD;          
	muD  = steady_muD;          
	UCD  = steady_UCD;          
	
	f_fcn   = 0;
	g_fcn   = 0;
	f_fcn_prime = 0;
	g_fcn_prime = 0;	

end;

%===================================================================%
%%%%            RUN                                              %%%%
%===================================================================%

shocks;
var epsilon_chi;
stderr 1;
end;

% steady(solve_algo  = 0, maxit = 10000);
steady(solve_algo  = 1, maxit = 50000);


check;

% Set seed for simulation
set_dynare_seed(092677);

% Produce simulation
% stoch_simul(order=1,periods=10000,nograph); 

% stoch_simul(order=1,periods=600, nograph, irf=10, loglinear);
