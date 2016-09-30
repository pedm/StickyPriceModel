function F = fbnd(x)
load TEMP.mat

alpha = 0.33;

% F(1) = (x(1)+1)*(10-x(1))*(1+x(2)^2)/(1+x(2)^2+x(2));
% F(2) = (x(2)+2)*(20-x(2))*(1+x(1)^2)/(1+x(1)^2+x(1));


g_fcn = 0;
g_fcn_prime = 0;
f_fcn = 0;
f_fcn_prime = 0;


%% Guess Two SS Values
g = x(1);
lambda = x(2);

%% Sticky Price Model (September Update)
ss_given_g_and_lambda;

%% Residuals: Equations 318 and 322
F(1) = (-Q + Lambda * ((g* (vartheta - 1) *YDW * alpha)/(mkup * KD * vartheta) + Q * (1 - delta)));
F(2) = (( chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) ) - ( (1/mkup) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L)));

