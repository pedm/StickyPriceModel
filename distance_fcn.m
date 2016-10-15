function [ f ] = distance_fcn( params_unbounded )

    % get Dynare structures
    global oo_ M_ options_ pvarcoirfs_clean
    % This ensures that we don't use the dr or irfs from the last stoch_simul()
    load level0workspace oo_ options_

    %% 1. Update parameters
    [ params ] = bounds( params_unbounded );
    set_param_value('eta', params(1) );
    set_param_value('gamma', params(2) );
    set_param_value('phi', params(3) );
    set_param_value('lambda_bar', params(4) );
    set_param_value('psi_N', params(5) );
    set_param_value('rhozeta', params(6) );
    set_param_value('rhozeta2', params(7) );
    set_param_value('sigmazeta', params(8) );
    % set_param_value('zetabar', params(9) ); % removed because it is set by ss solver
    set_param_value('rho_lambda', params(9) ); 

    try
        %% 2. Solve Again
        % Here we want to run the following dynare command:
        % stoch_simul(order=1,periods=600, irf=10, nograph, nodisplay, nocorr, nofunctions, nomoments, noprint, loglinear);
        % Dynare commands do not work in m files. Here's what it is in matlab:
        
        options_.irf = 11;
        options_.loglinear = 1;
        options_.nocorr = 1;
        options_.nodisplay = 1;
        options_.nograph = 1;
        options_.nomoments = 1;
        options_.noprint = 1;          % this turns off all display
        options_.order = 1;
        options_.periods = 600;
        var_list_=[];
        info = stoch_simul(var_list_); 
        
        %% 3. Compute Distance from IRFs
        % TODO: perhaps I can refactor this code so it's faster (for instance,
        % dont load the VAR irfs every single time. and dont compute so many
        % extra objects
        post_processing_irfs;                                                       % Create IRFs with trend
        post_processing_irfs_distance;                                              % Compute distance between model and VAR IRFs
        disp('ss and irfs found')
        
        if sum(oo_.steady_state == Inf) > 0
            % Dynare got a steady state with Inf values
            disp('Dynare got a steady state with Inf values. Why?')
            f = 10000000000;
        else
            f = irf_distance;
        end
    catch
        % disp('Error: No ss or no irfs found')
        
%         %% Add penalty
%         % I notice that csminwel often fails because its gradient is zero
%         % when it cannot find a steady state. Hopefully this creates a
%         % nonzero gradient
%         
%         compute_residuals;
%         
%         % Hmm. I'm not sure what to do with NaN
%         % should I replace with 0 or with something large?
%         how_many_nan = sum(isnan(z));
%         z(isnan(z)) = 1;
%         
%         sum_residuals = sum(abs(z));
        sum_residuals = 0;
        
        % Dynare threw an error. Apply large penalty
        f = 10000000000 + 10000*sum_residuals;
    end
end

