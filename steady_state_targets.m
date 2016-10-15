YD_ss = oo_.steady_state(strmatch('YDW' , M_.endo_names));
ND_ss = oo_.steady_state(strmatch('ND' , M_.endo_names));
CD_ss = oo_.steady_state(strmatch('CD' , M_.endo_names));
lambda_ss = oo_.steady_state(strmatch('lambda' , M_.endo_names));

disp(sprintf('ND/YD = %0.5g (target = 0.014)', ND_ss / YD_ss))
disp(sprintf('CD/YD = %0.5g (target = 0.75)', CD_ss / YD_ss))
disp(sprintf('lambda_ss = %0.5g', lambda_ss))



