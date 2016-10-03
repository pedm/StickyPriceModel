function [ f ] = fittingfunction10( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)

% WARNING: stoch_simul() will inherit whatever options were previously run
% PERHAPS: hard code these options_ in this function

    % get Dynare structures
    global oo_ M_ options_ pvarcoirfs_clean

    %% Set parameters, solve again, and plot
    
    % some of these restrictions arent entirely needed, but perhaps good
    % TODO: include gamma, phi, eta, lambda, psi_N (0 to Inf), rhozeta (0 to 1), sigmazeta (0 to inf), zetabar
    set_param_value('eta', p1 );
    set_param_value('gamma', p2 );
    set_param_value('phi', p3 );
    set_param_value('lambda_bar', p4 );
    set_param_value('psi_N', p5 );
    set_param_value('rhozeta', p6 );
    set_param_value('rhozeta2', p7 );
    set_param_value('sigmazeta', p8 );
    % set_param_value('zetabar', p9 );
    set_param_value('rho_lambda', p10 ); % woops. before I was setting rho_lambda = zetabar

    try
        % Dynare command - does not work in m file
        % stoch_simul(order=1,periods=600, irf=10, nograph, nodisplay, nocorr, nofunctions, nomoments, noprint, loglinear);
        % Here's what that dynare command translates to (assuming options_ already
        % defined)
        var_list_=[];
        info = stoch_simul(var_list_);
        
        % TODO: perhaps I can refactor this code so it's faster (for instance,
        % dont load the VAR irfs every single time. and dont compute so many
        % extra objects
        post_processing_irfs;                                                       % Create IRFs with trend
        post_processing_irfs_distance;                                              % Compute distance between model and VAR IRFs

        % TODO: eventually ill want to use the quadratic deviation of the
        % simulated irfs from the target irfs. perhaps see the code in Born and Pfeifer (2014)
        
        if sum(oo_.steady_state == Inf) > 0
            % Dynare got a steady state with Inf values
            disp('Dynare got a steady state with Inf values. Why?')
            % params
            % keyboard
            f = 10000000000;
        else
            f = irf_distance;
        end
    catch
        % params
        disp('Error: negative steady state')
        
        % Dynare threw a command. Apply large penalty
        f = 10000000000;
    end


    
end


