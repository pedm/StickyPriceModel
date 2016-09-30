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
    % Search over g, find other steady state values, then look at the residual of the final equation
    
    % 1. Function: residual given guess for growth rate g
    function f = subfunction(xx)

        g_fcn = 0;
        g_fcn_prime = 0;
        f_fcn = 0;
        f_fcn_prime = 0;
        
        % for i=1:cols
            
            %% Guess Two SS Values
            g = mintomod_ab(xx(1,1), phi*1.001, 2);
            lambda = mintomod_ab(xx(2,1), 0.0001, 0.9999);

            %% Sticky Price Model (September Update)
            ss_given_g_and_lambda;
            
            %% Residuals: Equations 318 and 322
            f2 = (( chi * GammaD * L^epsilon * (1/UCD) * (CD - GammaD * ( chi / (1+epsilon)) * L^(1+epsilon))^(-rho) ) - ( (1/mkup) * ((vartheta - 1)/vartheta) * (1 - alpha) * (YD/L)));
            f1 = (-Q + Lambda * ((g* (vartheta - 1) *YDW * alpha)/(mkup * KD * vartheta) + Q * (1 - delta)));
            
            
            % We get imaginary residuals if g, lambda, Z, K, V, or others are negative
            % To prevent negative VD, make phi small
%             if imag(f2) ~= 0
%                 disp('Imaginary residual')
%                 % keyboard
%                 
%                 % Not sure if this will work.... might still cause issues
%                 % because f is positive. I think that messes with root
%                 % finding algorithms (think: mean value theorem)
%                 f = [1000000; 1000000];
%             else
                 f = [f1; f2];
%             end
        % end
        
    end

    %% 2. Find growth rate g such that residual is zero
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
    
    opts = optimoptions(@lsqnonlin,'Display', 'iter-detailed');
    [x1b, RESNORM, ~, ~] = lsqnonlin(@fbnd,x0,lb, [Inf, 1], opts);    
%     RESNORM
%     F = fbnd(x1b)
    
    % TODO: if resnorm is large, try the algorithm again using the lower x0
    
    
  
    
% Other Solution Suggested by Matlab:::::::::::::::::::::::::::::::::::::::
% Set Equations and Inequalities as fmincon Constraints
% You can reformulate the problem and use fmincon as follows:
% 
% Give a constant objective function, such as @(x)0, which evaluates to 0 for each x.
% Set the fsolve objective function as the nonlinear equality constraints in fmincon.
% Give any other constraints in the usual fmincon syntax.
% For this example, write a function file for the nonlinear inequality constraint.
% 
% function [c,ceq] = fminconstr(x)
% 
% c = []; % no nonlinear inequality
% ceq = fbnd(x); % the fsolve objective is fmincon constraints
% Save this code as the file fminconstr.m on your MATLAB path.
% 
% Solve the constrained problem.
% 
% lb = [0,0]; % lower bound constraint
% rng default % reproducible initial point
% x0 = 100*randn(2,1);
% opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
% x = fmincon(@(x)0,x0,[],[],[],[],lb,[],@fminconstr,opts)

    
    %[xVal, fVal] = fzero(@subfunction, [phi*1.0001 2])
    
    % perhaps this fails because of complex f
    % [g_sol, rc]=csolve(@subfunction, g0, [], 1e-10, 1000);
    % This tells me whether it worked
    % rc
    
%     if rc == 4
%         disp('rc is 4 in steady state solver (maxit reached)')
%     end
    
    %g = g_sol(1,1)
    %lambda = mintomod_ab(g_sol(2,1), 0.01, 0.99)

     %           g = mintomod_ab(g_sol(1,1), phi*1.001, 2)
     %      lambda = mintomod_ab(g_sol(2,1), 0.0001, 0.9999);

     % g = x1b(1)
     g = 1.0118;
     zetabar = x1b(1)
     lambda = x1b(2)
     
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

