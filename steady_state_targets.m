if options_.logged_steady_state == 1
    sstate = exp(oo_.steady_state);
else
    sstate = oo_.steady_state;
end

YD_ss = sstate(strmatch('YDW' , M_.endo_names));
ND_ss = sstate(strmatch('ND' , M_.endo_names));
CD_ss = sstate(strmatch('CD' , M_.endo_names));
lambda_ss = sstate(strmatch('lambda' , M_.endo_names));

disp(sprintf('\nSteady State Targets:'));
disp(sprintf('ND/YD = %0.5g (target = 0.014)', ND_ss / YD_ss))
disp(sprintf('CD/YD = %0.5g (target = 0.75)', CD_ss / YD_ss))
disp(sprintf('lambda_ss = %0.5g', lambda_ss))


zeta_bar_value = param_value( 'zeta_bar' );
lambda_bar_value = param_value( 'lambda_bar' );
chi_value = param_value( 'chi' );

disp(sprintf('\nParameters set by SS Solver:'));
disp(sprintf('zeta_bar = %0.5g', zeta_bar_value ))
disp(sprintf('lambda_bar = %0.5g', lambda_bar_value ))
disp(sprintf('chi = %0.5g', chi_value ))



