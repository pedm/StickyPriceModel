function [ f ] = distance_fcn( params_unbounded )

% WARNING: stoch_simul() will inherit whatever options were previously run
% PERHAPS: hard code these options_ in this function

    % get Dynare structures
    global oo_ M_ options_ 

    %% Change parameters, solve again, and plot
    [ params ] = bounds( params_unbounded );
    params
    set_param_value('eta', params(1) );
    set_param_value('gamma', params(2) );
    set_param_value('phi', params(3) );
    
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
    f = irf_distance;
    
end

