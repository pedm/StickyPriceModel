function [ ss ] = ss_solver_flex_price_model( param_names, params)
    
    %% Load parameters being used in mod file
    
    % Assign variables (needed for dynamic assignment using eval)
    beta    = [];
    alpha   = [];
    epsilon = [];
    rho     = [];
    varphi = [];
    delta   = [];
    chi     = [];
    vartheta = [];
    gamma  = [];
    phi    = [];
    eta    = [];
    lambda = [];
    LS = [];
    rhozeta = [];
    M       = [];
    psi_N   = [];
    psi_I   = [];
    gg      = [];

    % Load variable parameters
    for ii = 1:length(params)
        cmd = sprintf([param_names(ii,:), ' = %0.6f;'], params(ii));
        eval(cmd);
    end
    clear gg

    %% Solve for Steady State g
    % Search over g, find other steady state values, then look at the residual of the final equation
    
    function f = subfunction(xx)

        [~,cols]=size(xx);
        
        g_fcn = 0;
        g_fcn_prime = 0;
        f_fcn = 0;
        f_fcn_prime = 0;
        
        for i=1:cols
            
            g =  xx(1,i);

            %% Flex Price Model (August Update)
            
            % Use eqns 1 and 2
            VD = ( (g-phi) * (g - phi + phi*lambda) ) / ( phi * lambda * (g - phi*lambda) );
            ZD = VD * phi * (1-lambda) / (g-phi);
            
            % Eqn 19
            zeta = 1;
            
            % Eqn 3
            ND = (      VD / (zeta * ZD ^ (1-eta))      )^(1/eta);
            
            % Eqn 11
            Lambda = beta * g^(-rho);
            
            % Eqn 7
            J =  1 / ( eta * zeta * ( ZD / ND )^(1-eta)  );

            % Eqn 4
            H = (J -  ((1 - lambda) * phi * Lambda * J ) ) * (1/lambda);
            
            % Eqn 5
            Pi = H - ( phi * Lambda * H ) ;
            
            % Eqn 6
            YDW = Pi / ( (1/vartheta) * (1/M) );
            
            % Eqn 8
            YD = YDW;
            
            % Eqn 17
            Q = 1;
            
            % Eqn 18
            block = (g* (vartheta - 1) *YDW * alpha)/(M * vartheta);
            KD = Lambda * block / ( Q - Lambda * Q * (1 - delta) ); 
            
            % Eqn 9
            L = (  YDW / (    (KD/g) ^ alpha    )   )^(1/(1 - alpha));
            
            % Eqn 16
            ID = KD - (1-delta) * (KD/ g);
            
            % Eqn 10
            % CD = YD - (1 + g_fcn) * ID + ND;
            CD = YD - (1 + g_fcn) * ID - ND;
            
            % Eqn 15
            GammaD = CD * g^((gamma - 1)/gamma);
            
            % Eqn 13
            % muD = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) - ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
            AA = beta * (1-gamma) * ( (g^(-rho)) * (CD * g/ GammaD)^gamma );
            BB = - ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
            % muD = AA*muD + BB;
            muD = BB / (1-AA);
            
            % Eqn 12
            UCD = ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) + muD * gamma * (GammaD / ( CD * g )) ^ (1-gamma);
            
            % Eqn 14
            f = ( chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) ) - ( (1/M) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L));
       
        end
        
    end


g0 = 1.0001;
[g_sol, rc]=csolve(@subfunction, g0, [], 1e-12, 800);rc
g = g_sol;

% solve again (same as above)
% Use eqns 1 and 2
VD = ( (g-phi) * (g - phi + phi*lambda) ) / ( phi * lambda * (g - phi*lambda) );
ZD = VD * phi * (1-lambda) / (g-phi);

% Eqn 19
zeta = 1;

% Eqn 3
ND = (      VD / (zeta * ZD ^ (1-eta))      )^(1/eta);

% Eqn 11
Lambda = beta * g^(-rho);

% Eqn 7
J =  1 / ( eta * zeta * ( ZD / ND )^(1-eta)  );

% Eqn 4
H = (J -  (1 - lambda) * phi * Lambda * J ) / lambda;

% Eqn 5
Pi = H - phi * ( Lambda * H ) ;

% Eqn 6
YDW = Pi / ( (1/vartheta) * (1/M) );

% Eqn 8
YD = YDW;

% Eqn 17
Q = 1;

% Eqn 18
block = (g* (vartheta - 1) *YDW * alpha)/(M * vartheta);
KD = Lambda * block / ( Q - Lambda * Q * (1 - delta) ); 

% Eqn 9
L = (  YDW / (    (KD/g) ^ alpha    )   )^(1/(1 - alpha));

% Eqn 16
ID = KD - (1-delta) * (KD/ g);

% Eqn 10
CD = YD - (1 + g_fcn) * ID - ND;

% Eqn 15
GammaD = CD * g^((gamma - 1)/gamma);

% Eqn 13
% muD = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) - ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
AA = beta * (1-gamma) * ( g^(-rho) * (CD * g/ GammaD)^gamma );
BB = - ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
% muD = AA*muD + BB;
muD = BB / (1-AA);

% Eqn 12
UCD = ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) + muD * gamma * (GammaD / ( CD * g )) ^ (1-gamma);


ss.g          = g;
ss.VD          = VD ;
ss.ZD          = ZD ;
ss.zeta          = zeta ;
ss.ND          = ND ;
ss.Lambda          = Lambda ;
ss.J          = J ;
ss.H          = H ;
ss.Pi          = Pi ;
ss.YDW          = YDW ;
ss.YD          = YD ;
ss.Q          = Q ;
ss.KD          = KD ;
ss.L          = L ;
ss.ID          = ID ;
ss.CD          = CD ;
ss.GammaD          = GammaD ;
ss.muD          = muD ;
ss.UCD          = UCD ;
    
ss

% collect into structure
steady_g          = g;
steady_VD          = VD ;
steady_ZD          = ZD ;
steady_zeta          = zeta ;
steady_ND          = ND ;
steady_Lambda          = Lambda ;
steady_J          = J ;
steady_H          = H ;
steady_Pi          = Pi ;
steady_YDW          = YDW ;
steady_YD          = YD ;
steady_Q          = Q ;
steady_KD          = KD ;
steady_L          = L ;
steady_ID          = ID ;
steady_CD          = CD ;
steady_GammaD          = GammaD ;
steady_muD          = muD ;
steady_UCD          = UCD ;
    
% Save it. Load into dynare later
save ss_flexprice.mat steady_*

end

