clear;clc;


% PARAMS

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

% STEADY STATE

% 
%          g =  1.2453;
%         VD =  1.4459;
%         ZD =  3.3913;
%       zeta =  1;
%         ND =  0.2561;
%     Lambda =  0.7709;
%          J =  0.5368;
%          H =  2.0161;
%         Pi =  0.6174;
%        YDW =  2.0247;
%         YD =  2.0247;
%          Q =  1;
%         KD =  0.9533;
%          L =  3.2690;
%         ID =  0.2644;
%         CD =  2.0165;
%     GammaD =  0.4035;
%        muD =  -3.4178;
%        UCD =  0.5200;
       

         g = 1.2733
        VD = 1.6241
        ZD = 3.5239
      zeta = 1
        ND = 0.3370
    Lambda = 0.7539
         J = 0.6288
         H = 2.4478
        Pi = 0.7869
       YDW = 2.5806
        YD = 2.5806
         Q = 1
        KD = 1.1574
         L = 4.3144
        ID = 0.3393
        CD = 1.9043
    GammaD = 0.3237
       muD = -5.8481
       UCD = 0.5601


f_fcn   = 0;
g_fcn   = 0;
f_fcn_prime = 0;
g_fcn_prime = 0;


% shock
epsilon_chi=0;

% set gg based on steady state
gg      = g;

%=========================ENDOGENOUS VARIABLES============================%

% 1. Evolution of number of firm varieties
check.eq1 = g - phi* (1 + lambda * (ZD + VD - 1)) 

% 2. Evolution of technological fronteir Z
check.eq2 = ZD * g - phi * ( ZD + (1-lambda) * VD ) 

% 3. Innovators' production function. Q: Why is this not the supply curve of new firms?
% PAT'S MODIFICATION (I replace exp(zeta) with zeta to account for zeta being 1 in steady state)
check.eq3 = VD - zeta * ZD ^ (1-eta) * ND ^ eta 

% PREVIOUS: 3. Supply curve of new firms
% J = ((1/zeta)^(1 / eta) * (1/eta) *(1 /LS )^((1 - eta)/eta)*(VD/ZD)^((1 - eta)/eta)) 

% 4. Euler equation for entrepreneurs
check.eq4 = J -  (lambda * H + (1 - lambda) * phi * Lambda * J );

% 5. Price of adopted technology H
check.eq5 = H - (Pi + phi * ( Lambda * H ) );

% 6. Profits per period
check.eq6 = Pi - (1/vartheta) * (1/M) * YDW 

% 7. Innovators' FOC wrt N
% PAT'S MODIFICATION (I replace exp(zeta) with zeta to account for zeta being 1 in steady state)
check.eq7 = eta * J * zeta * ( ZD / ND )^(1-eta) - ( 1 + f_fcn_prime *  (ND * g /  ND) +  f_fcn - Lambda * f_fcn_prime * (ND * g / ND )^2 );

% PREVIOUS: 7. Final output used by the entrepreneurial sector
% ND = eta * VD * J 

% 8. Aggregate production function
check.eq8 = YD - YDW 

% 9. Aggregate production function W
check.eq9 = YDW - ((KD / g)^alpha) * L^(1-alpha) 

% 10. Resource constraint, with adjustment cost
check.eq10 = YD - (CD + (1 + g_fcn) * ID + ND )

% 11. HH's stochastic discount factor
check.eq11 = Lambda - ((beta * UCD) / UCD) * g^(-rho) 

% 12. Marginal utility of consumption
check.eq12 = UCD -( ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) + muD * gamma * (GammaD / ( CD * g )) ^ (1-gamma) );

% 13. Lagrange multiplier on labor disutility law of motion (new) (equation 258)
% Albert: equation 258 needs to be corrected - it will need a “g” term as well.
% To use dynare notation, here's the change C / Gamma = CD * g / GammaD
check.eq13 = muD - ( beta * (1-gamma) * ( g^(-rho) * muD * (CD * g/ GammaD)^gamma ) - ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon) );

% 14. Labor market check.equilibrium (eqn 259)
check.eq14 = chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) - (1/M) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L) 

% 15. Labor disutility term
check.eq15 = GammaD - (CD^gamma) * (GammaD / g )^(1-gamma) 

% 16. Capital evolution
check.eq16 = KD - ( (1-delta) * (KD / g) + ID  );

% 17. Q-equation (capital producers)
check.eq17 = Q - ( 1 + g_fcn + ((ID * g) / ID) * g_fcn_prime - Lambda * ((ID * g) / ID)^2 * g_fcn_prime );

% PREVIOUS: 17. Q-equation
% Q = ((ID * g * (I_K_ss^-1)) / (KD )  )^varphi 

% 18. equation 263
check.eq18 = Q - Lambda * ((g* (vartheta - 1) *YDW * alpha)/(M * KD * vartheta) + Q * (1 - delta)) 

% 19. Exogenous shock to entrepreneurs' production function
check.eq19 = log(zeta) - rhozeta * log(zeta) + 0.1 * epsilon_chi 

% 20. Aggregate R&D expenditure
% check.eq20 = RD - J * VD 

% ADJUSTMENT COST FUNCTIONS
check.eq21 = f_fcn - (psi_N / 2) * ((ND * g / ND ) - gg)^2 
check.eq22 = f_fcn_prime - psi_N * ((ND * g / ND ) - gg) 

check.eq23 = g_fcn - (psi_I / 2) * ((ID * g / ID ) - gg)^2 
check.eq24 = g_fcn_prime - psi_I * ((ID * g / ID ) - gg)


check

% Issues:
% check.eqn 4, 5 (big), 10, 12, 13 (biggest), 16, 

