g = 0.4
lambda = 0.1

ss_given_g_and_lambda;
f2 = (( chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) ) - ( (1/mkup) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L)))
f1 = (-Q + Lambda * ((g* (vartheta - 1) *YDW * alpha)/(mkup * KD * vartheta) + Q * (1 - delta)))


g = 10
lambda = 0.1

ss_given_g_and_lambda;
f2 = (( chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) ) - ( (1/mkup) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L)))
f1 = (-Q + Lambda * ((g* (vartheta - 1) *YDW * alpha)/(mkup * KD * vartheta) + Q * (1 - delta)))

g = 0.4
lambda = 0.9

ss_given_g_and_lambda;
f2 = (( chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) ) - ( (1/mkup) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L)))
f1 = (-Q + Lambda * ((g* (vartheta - 1) *YDW * alpha)/(mkup * KD * vartheta) + Q * (1 - delta)))



jj = 1;
FF = [];
GG = [];
LAMBDA = [];
for g = 0.9:.05:1.3
    iii = 1;
    for lambda = .1:.1:0.9
        ss_given_g_and_lambda;
        f2 = (( chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) ) - ( (1/mkup) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L)));
        f1 = (-Q + Lambda * ((g* (vartheta - 1) *YDW * alpha)/(mkup * KD * vartheta) + Q * (1 - delta)));
        
        GG(iii, jj) = g;
        LAMBDA(iii, jj) = lambda;
        FF(iii, jj) = abs(f1) + abs(f2);
        iii = iii+1;
    end
    jj = jj+1;
end

surf(GG,LAMBDA,FF)

