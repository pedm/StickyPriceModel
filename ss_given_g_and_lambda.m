%% Sticky Price Model (September Update)
% Given g and lambda, solve for remaining variables in steady state

%% Obvious SS Solutions
Lambda = beta * g^(-rho);
mkup = mkup_ss;
zeta = 1;
Q = 1;

%% Other Variables
% Equation 303
ZD = 1 + (g - phi) / lambda*phi;            

% Equation 304
VD = g*ZD - phi*ZD;

% Equation 311
M = (lambda / lambda_bar) ^ (1/rho_lambda);

% Equations 310 and 306
% Solve jointly for J and H
J = (1/(1 - phi*Lambda)) * (-M + (lambda * M^(1 - rho_lambda)) / (rho_lambda * lambda_bar));

H = J + (M^(1 - rho_lambda)) / (rho_lambda * lambda_bar * phi * Lambda);

% Equation 307
Pi = H - phi * Lambda * H;

% Equation 308
YDW = vartheta * mkup * Pi;
YD = YDW;

% Equations 305 and 309
% Solve jointly for N and K
ND = Lambda * J * VD;
KD = (zetabar * zeta * (ZD/VD)) ^ (1/eta) * ND * g;

% Equation 312
L = ( YDW * (g/KD) ^ alpha  )^(1/(1 - alpha));

% Equation 320
ID = KD - ((1-delta)/g)*KD;

% Equation 314
CD = YD - ID - ND - (ZD-1)*M;

% Equation 319
GammaD = CD * g^((gamma - 1)/gamma);

% Equation 317
% Original:
% muD = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) + ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
% To make muD positive, replace - with +
% muD = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) - ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
AA = beta * (1-gamma) * ( (g^(-rho)) * (CD * g/ GammaD)^gamma );
BB = ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
% muD = AA*muD + BB;
muD = BB / (1-AA);

% Equation 316
% Original:
% UCD = ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) + muD * gamma * (GammaD(-1) / ( CD * g(-1) )) ^ (1-gamma);
% To make muD positive, multiply by -1
UCD = ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) - muD * gamma * (GammaD / ( CD * g )) ^ (1-gamma);
