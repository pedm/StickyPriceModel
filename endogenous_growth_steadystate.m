function[ys,check]=endogenous_growth_steadystate(ys,exe)
    % Steady state file for the Flex Price Model of Endogenous Growth
    % August/Sept 2016 Update
    % Includes Jaimovich-Rebelo Preferences, a simplified innovation sector, and adjustment costs
    % Patrick Moran (based on example code by Steffen Esser)
    
    global M_
    
    %% Name the structural parameters 
    % Matlab requires each variable to be named somewhere in the function. 
    % Otherwise it throws an error when using eval() to define the variables.
    % This has to do with static vs dynamic assignment. Super boring.
    
    beta      = [];
    alpha     = [];
    epsilon   = [];
    rho       = [];
    varphi    = [];
    delta     = [];
    chi       = [];
    vartheta  = [];
    gamma     = [];
    phi       = [];
    eta       = [];
    lambda    = [];
    rhozeta   = [];
    M         = [];
    psi_N     = [];
    psi_I     = [];
    gg        = [];

    %% Load the values of the structural parameters in a loop.
    
    NumberOfParameters = M_.param_nbr;                            % Number of structural parameters.
    for i = 1:NumberOfParameters                                  % Loop...
      paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
      eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
    end                                                           % End of the loop.  
    % TODO what does this do?
    check = 0;
    
    %% Solve for Steady State:
    % Search over g, find other steady state values, then look at the residual of the final equation
    
    % 1. Function: residual given guess for growth rate g
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

    %% 2. Find growth rate g such that residual is zero
    g0 = 1.0001;
    [g_sol, rc]=csolve(@subfunction, g0, [], 1e-12, 800);
    % This tells me whether it worked
    % rc
    g = g_sol;

    %% 3. Given solution g, find the remaining steady state variables (same equations as above)
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
    % muD = beta * (1-gamma) * ( g^(-rho) * muD(+1) * (CD(+1) * g/ GammaD)^gamma ) + ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
    AA = beta * (1-gamma) * ( g^(-rho) * (CD * g/ GammaD)^gamma );
    BB = ((CD - GammaD*( chi / (1+epsilon)) * L^(1+epsilon))^(-rho)) * ( chi / (1+epsilon)) * L^(1+epsilon);
    % muD = AA*muD + BB;
    muD = BB / (1-AA);

    % Eqn 12
    UCD = ( CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon) ) ^ (-rho) + -1*muD * gamma * (GammaD / ( CD * g )) ^ (1-gamma);
    
    %% State the remaining steady state variables

    % XD =  ( Lambda * g * ( J * VD + XD ) );
    XD =  ( Lambda * g * ( J * VD )  ) / (1 - Lambda * g ) ;
    SD = Q * KD + H  +  J * ( ZD + VD - 1 )  + XD;
    RD = J * VD;
    
    f_fcn   = 1;
	g_fcn   = 1;
	f_fcn_prime = 1;
	g_fcn_prime = 1;	
    
    %% Save local parameter values to the M_ global
    % This is useful for any parameters that depend on the steady state
    
    % Currently, only parameter that depends on the ss is gg
    gg = g_sol;
    
    for iter = 1:length(M_.params) %update parameters set in the file
      eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
    end

    %% Output the values of the endogenous vars at steady state

    NumberOfEndogenousVariables = M_.endo_nbr;                    % Number of endogenous variables.
    ys = zeros(NumberOfEndogenousVariables,1);                    % Initialization of ys (steady state).
    for i = 1:NumberOfEndogenousVariables                         % Loop...
      varname = deblank(M_.endo_names(i,:));                      %    Get the name of endogenous variable i.                     
      eval(['ys(' int2str(i) ') = ' varname ';']);
      % ys(i)=0;%    Or just save every ss as zero
    end   
    
end

