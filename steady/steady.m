function steady = steady(param);

% find steady state given params

beta = param.beta;
alpha= param.alpha;
epsilon = param.epsilon;
rho = param.rho;
chi = param.chi;

vartheta = param.vartheta;
delta    = param.delta;

gamma = param.gamma;

phi    = param.phi;
eta    = param.eta;
LS     = param.LS;
lambda = param.lambda; 
M      = param.M;

% now solve steady state; search over g.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOLVE STEADY STATE  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

save ss_solution.mat g VD  ZD  zeta  ND  Lambda  J  H  Pi  YDW  YD  Q  KD  L  ID  CD  GammaD  muD  UCD

% collect into structure
steady.g          = g;
steady.VD          = VD ;
steady.ZD          = ZD ;
steady.zeta          = zeta ;
steady.ND          = ND ;
steady.Lambda          = Lambda ;
steady.J          = J ;
steady.H          = H ;
steady.Pi          = Pi ;
steady.YDW          = YDW ;
steady.YD          = YD ;
steady.Q          = Q ;
steady.KD          = KD ;
steady.L          = L ;
steady.ID          = ID ;
steady.CD          = CD ;
steady.GammaD          = GammaD ;
steady.muD          = muD ;
steady.UCD          = UCD ;


% % % collect in X, then names
% % X(1,1) = g;
% % X(2,1) = ZD;
% % X(3,1) = VD;
% % X(4,1) = J;
% % X(5,1) = H;
% % X(6,1) = Lambda;
% % X(7,1) = L;
% % X(8,1) = KD;
% % X(9,1) = YD;
% % X(10,1)= CD;
% % X(11,1)= ID;
% % X(12,1)= ND;
% % X(13,1)= CD/YD;
% % X(14,1)= ID/YD;
% % X(15,1)= ND/YD;
% % X(16,1)= g^4-1;
% % 
% % names = {'g','ZD','VD','J','H','Lambda','L','KD','YD','CD','ID','ND','C/Y','I/Y','N/Y','annual g'};


end





