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
    rhozeta2  = [];
    sigmazeta = [];
    zetabar   = [];
    mkup_ss   = [];
    psi_N     = [];
    psi_I     = [];
    gg        = [];
    R_nom_ss  = [];
    gamma_pi  = [];
    gamma_y  = [];
    omega  = [];
    theta  = [];
    lambda_bar = [];
    rho_lambda = [];

    % Name the endogenous variables
    Lambda = [];
    mkup = [];
    zeta = [];
    Q = [];
    ZD = [];
    VD = [];
    M = [];
    J = [];
    H = [];
    Pi = [];
    YDW = [];
    YD = [];
    ND = [];
    KD = [];
    L = [];
    ID = [];
    CD = [];
    GammaD = [];
    AA = [];
    BB = [];
    muD = [];
    UCD = [];

    %% Load the values of the structural parameters in a loop.
    
    NumberOfParameters = M_.param_nbr;                            % Number of structural parameters.
    for i = 1:NumberOfParameters                                  % Loop...
      paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
      eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
    end                                                           % End of the loop.  
    % TODO what does this do?
    check = 0;
    

    
    %% Solve for Steady State:
    % Search over lambda and zetabar, find other steady state values, then look at the residual of the final equation
    
    %% 1. Find lambda and zetabarsuch that residual is zero
    % g0 = [1.01; 0];
    
    save TEMP.mat
    % Usually the best choice of bounds
    x0 = [1.01, 0.5];

    % Sometimes this works better
    % x0 = [1.01, 0.1];
    % But it usually results in lambda in ss being VERY close to zero
    
%     opts = optimoptions(@fsolve,'Display', 'iter','Algorithm','trust-region-dogleg');
%     [x1, ff] = fsolve(@fbnd,x0, opts)

    % Solve for the roots of fbnd, while also imposing a lower bound
    % Point of the lower bound is to avoid complex results
    % For more info, look in matlab documentation on "Nonlinear Systems
    % with Constraints"
    
    % lb when guessing g
    % lb = [phi,0];
    
    % lb when guessing zetabar
    lb = [0.01, 0];
    
    opts = optimoptions(@lsqnonlin,'Display', 'iter-detailed'); % debugging
    opts = optimoptions(@lsqnonlin,'Display', 'off'); % estimation

    % rng(23942374) % TODO: or is this why the results changed slightly?
    
    [x1b, RESNORM, ~, ~] = lsqnonlin(@fbnd,x0,lb, [Inf, 1], opts);    
    % RESNORM
    % F = fbnd(x1b)
    
    % if resnorm is large, try the algorithm using a lower starting value
    if RESNORM > 1e-10
        % disp('resnorm too high. try low x0')
        try
            x0 = [0.4, 0.1];
            [x1b_2, RESNORM_2, ~, ~] = lsqnonlin(@fbnd,x0,lb, [Inf, 1], opts);
            if RESNORM_2 < RESNORM
                x1b = x1b_2;
                RESNORM = RESNORM_2;
            end
        end
    end
    
    % if resnorm is still too large, try the algorithm using a high starting value
    if RESNORM > 1e-10
        % disp('resnorm too high. try high x0')
        try
            x0 = [20, 0.9];
            [x1b_2, RESNORM_2, ~, ~] = lsqnonlin(@fbnd,x0,lb, [Inf, 1], opts);
            if RESNORM_2 < RESNORM
                x1b = x1b_2;
                RESNORM = RESNORM_2;
            end
        end
    end

    

     g = 1.0118;
     zetabar = x1b(1);
     lambda = x1b(2);
     
    %% 3. Given solution, find the remaining steady state variables (same equations as above)
    ss_given_g_and_lambda;

    %% State the remaining steady state variables

    XD =  ( Lambda * g * ( J * VD )  ) / (1 - Lambda * g ) ;
    SD = Q * KD + H  +  J * ( ZD + VD - 1 )  + XD;
    RD = ND;
    
    f_fcn   = 1;
	g_fcn   = 1;
	f_fcn_prime = 1;
	g_fcn_prime = 1;	
    
    % Zero inflation steady state
    pi = 1;
    pi_star = 1;
    x1D = UCD*(1/mkup)*YD / (1 - theta*beta*g^(1-rho)*pi^omega);
    x2D = UCD*YD / (1 - theta*beta*g^(1-rho)*pi^(omega-1));
    R_nom = pi / Lambda;
    
    %% Save local parameter values to the M_ global
    % This is useful for any parameters that depend on the steady state
    
    % Currently, only parameter that depends on the ss is gg
    gg = g;
    R_nom_ss = R_nom;
    
    if gg>=1.5
        disp('g in steady state > 1.5')
    end
    
    for iter = 1:length(M_.params) %update parameters set in the file
        eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
    end
    
    %% Output the values of the endogenous vars at steady state

    NumberOfEndogenousVariables = M_.endo_nbr;                    % Number of endogenous variables.
    ys = zeros(NumberOfEndogenousVariables,1);                    % Initialization of ys (steady state).
    for i = 1:NumberOfEndogenousVariables                         % Loop...
      varname = deblank(M_.endo_names(i,:));                      %    Get the name of endogenous variable i.  
      
      % Hardcode ss = 1 for AUX_ENDO_LAG_... variable (appears if using AR(2))
      if length(varname) >= 8 && strcmp(varname(1:8), 'AUX_ENDO')
          eval(['ys(' int2str(i) ') = 1;']);
      else
          eval(['ys(' int2str(i) ') = ' varname ';']);
      end
      % ys(i)=0;%    Or just save every ss as zero
    end   
    
end

