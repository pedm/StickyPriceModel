function [ys,check] = endogenous_growth_sticky_steadystate(ys,exe);

% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%
% Output: 
%   - ys        [vector] vector of steady state values fpr the the endogenous variables
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impos restriction on parameters)
    
    
global M_ CD GammaD epsilon L rho beta gamma g mkup_ss vartheta alpha YD

    
alpha = [];
beta  = [];
gamma = []; % these 3 are needed so they don't clash with functions with the same name - in the future better to use different parameter names: bet, alph, gam

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end

% initialize indicator (check=failure indicator -- indicates ss computation
% failed)
check = 0;



%%%% ---------------------- %%%    
%% Solve for Steady State %%%%%
%%% ------------------------%%%
% solve for YD, KD
L     = L_ss;
g     = g_ss;
lambda= lambda_ss;
LAMBDA = beta*g^(-rho);
YD_KD_g= ( LAMBDA^(-1) + delta - 1 )/( (vartheta-1)/vartheta * 1/mkup_ss * alpha   ) ;
KD_g   = L * (YD_KD_g)^(-1/(1-alpha));
KD = KD_g * g;
YD = (KD/g)^alpha * L^(1-alpha);
PI = (1/vartheta)*(1/mkup_ss)*YD;
H  = (1/(1-phi*LAMBDA))*PI;
YDW = YD;

% solve for M and lambdabar (following the notes)
% M_coef_num   = 1-phi*( (1-lambda)*LAMBDA + lambda);
% M_coef_denom = 1-phi*( (1-lambda)*LAMBDA + lambda*rho_lambda);

denom = 1 - phi*LAMBDA*(1-lambda);
M_coef_num = 1 - phi * lambda * LAMBDA * (1/denom);
M_coef_denom = 1 - rho_lambda * lambda * phi * LAMBDA * (1/denom);

M = (M_coef_num/M_coef_denom) * rho_lambda * phi * LAMBDA * lambda * H;
J = (1/(1-phi*LAMBDA*(1-lambda)))*(-M + phi*LAMBDA*lambda*H);
lambda_bar = lambda/(M^rho_lambda);

ZD = 1 + (g-phi)/(lambda*phi);
VD = (g-phi)*ZD;
zeta_bar = VD^(1-eta) * (KD/g)^eta * (LAMBDA*J)^(-eta) * ZD^(-1) ;

ND = VD^(1/eta) * (zeta_bar*ZD)^(-1/eta) * KD/g;


% solve for the rest
ID = ( 1 - (1-delta)/g ) * KD;
CD = YD - ID - ND - (ZD-1)*M;
GammaD = CD * g^((gamma-1)/gamma);


% next solve for labor disutility parameter chi

% set starting value for chi
yy0(1,1) = 0.0001;

% solve for chi (see subfunction below)
[yy, rc]=csolve(@subfunction, yy0, [], 1e-12, 1000);rc
chi = yy
W  = ( CD - GammaD * (chi/(1+epsilon)) * L^(1+epsilon) )^(-rho);
kappa_mu = beta*(1-gamma)*g^(-rho)*(CD*g/GammaD)^gamma;
muD = (1/(1-kappa_mu))*W*(chi/(1+epsilon))*L^(1+epsilon);
UCD= W - muD*gamma*( GammaD/(CD*g) )^(1-gamma);

    

% sol. for adjustment cost functions
f_fcn   = 1;
g_fcn   = 1;
f_fcn_prime = 1;
g_fcn_prime = 1;	
    
% sol. for sticky-price variables (zero inflation steady state)
mkup = mkup_ss;
pi = 1;
pi_star = 1;
x1D = UCD*(1/mkup)*YD / (1 - theta*beta*g^(1-rho)*pi^omega);
x2D = UCD*YD / (1 - theta*beta*g^(1-rho)*pi^(omega-1));
R_nom = pi / LAMBDA;
    
Q = 1;
zeta = 1;

%%%%---------------------------%%%%
%% End Own Steady State Solution %%
%%%%---------------------------%%%%


for iter = 1:length(M_.params) %update parameters set in the file
    c = [ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ];
    % disp(c)
    eval(c)
end



NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end   
    
  
   
    
end



%%-------------------------%%%%
%%   subfunction of chi    %%%%
%%-------------------------%%%%
    function ff = subfunction(yy)
        global CD GammaD epsilon L rho beta gamma g mkup_ss vartheta alpha YD
        [~,cols]=size(yy);
         
        for i=1:cols            
            
            chi =  yy(1,i);   
            
            W  = ( CD - GammaD * (chi/(1+epsilon)) * L^(1+epsilon) )^(-rho);
            kappa_mu = beta*(1-gamma)*g^(-rho)*(CD*g/GammaD)^gamma;
            muD = -(1/(1-kappa_mu))*W*(chi/(1+epsilon))*L^(1+epsilon);
            UCD= W + muD*gamma*( GammaD/(CD*g) )^(1-gamma);
  
            ff(1,i) =  chi*W*GammaD*L^epsilon - UCD*(1/mkup_ss)*((vartheta-1)/vartheta)*(1-alpha)*YD/L;           
       
        end        
    end
%%-----------------------------%%%%
%%   end subfunction of chi    %%%%
%%-----------------------------%%%%









