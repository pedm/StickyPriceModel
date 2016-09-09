%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'endogenous_growth';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('endogenous_growth.log');
M_.exo_names = 'epsilon_chi';
M_.exo_names_tex = '{\epsilon}^{\chi}';
M_.exo_names_long = 'epsilon_chi';
M_.endo_names = 'g';
M_.endo_names_tex = '{g}';
M_.endo_names_long = 'g';
M_.endo_names = char(M_.endo_names, 'ZD');
M_.endo_names_tex = char(M_.endo_names_tex, '{Z_D}');
M_.endo_names_long = char(M_.endo_names_long, 'ZD');
M_.endo_names = char(M_.endo_names, 'VD');
M_.endo_names_tex = char(M_.endo_names_tex, '{V_D}');
M_.endo_names_long = char(M_.endo_names_long, 'VD');
M_.endo_names = char(M_.endo_names, 'J');
M_.endo_names_tex = char(M_.endo_names_tex, '{J}');
M_.endo_names_long = char(M_.endo_names_long, 'J');
M_.endo_names = char(M_.endo_names, 'H');
M_.endo_names_tex = char(M_.endo_names_tex, '{H}');
M_.endo_names_long = char(M_.endo_names_long, 'H');
M_.endo_names = char(M_.endo_names, 'Pi');
M_.endo_names_tex = char(M_.endo_names_tex, '{\Pi}');
M_.endo_names_long = char(M_.endo_names_long, 'Pi');
M_.endo_names = char(M_.endo_names, 'ND');
M_.endo_names_tex = char(M_.endo_names_tex, '{N_D}');
M_.endo_names_long = char(M_.endo_names_long, 'ND');
M_.endo_names = char(M_.endo_names, 'YDW');
M_.endo_names_tex = char(M_.endo_names_tex, '{Y^W_D}');
M_.endo_names_long = char(M_.endo_names_long, 'YDW');
M_.endo_names = char(M_.endo_names, 'YD');
M_.endo_names_tex = char(M_.endo_names_tex, '{Y_D}');
M_.endo_names_long = char(M_.endo_names_long, 'YD');
M_.endo_names = char(M_.endo_names, 'CD');
M_.endo_names_tex = char(M_.endo_names_tex, '{C_D}');
M_.endo_names_long = char(M_.endo_names_long, 'CD');
M_.endo_names = char(M_.endo_names, 'Lambda');
M_.endo_names_tex = char(M_.endo_names_tex, '{\Lambda}');
M_.endo_names_long = char(M_.endo_names_long, 'Lambda');
M_.endo_names = char(M_.endo_names, 'UCD');
M_.endo_names_tex = char(M_.endo_names_tex, '{U_{CD}}');
M_.endo_names_long = char(M_.endo_names_long, 'UCD');
M_.endo_names = char(M_.endo_names, 'muD');
M_.endo_names_tex = char(M_.endo_names_tex, '{{\mu}_{D}}');
M_.endo_names_long = char(M_.endo_names_long, 'muD');
M_.endo_names = char(M_.endo_names, 'L');
M_.endo_names_tex = char(M_.endo_names_tex, '{L}');
M_.endo_names_long = char(M_.endo_names_long, 'L');
M_.endo_names = char(M_.endo_names, 'GammaD');
M_.endo_names_tex = char(M_.endo_names_tex, '{\Gamma_D}');
M_.endo_names_long = char(M_.endo_names_long, 'GammaD');
M_.endo_names = char(M_.endo_names, 'KD');
M_.endo_names_tex = char(M_.endo_names_tex, '{K_D}');
M_.endo_names_long = char(M_.endo_names_long, 'KD');
M_.endo_names = char(M_.endo_names, 'Q');
M_.endo_names_tex = char(M_.endo_names_tex, '{Q}');
M_.endo_names_long = char(M_.endo_names_long, 'Q');
M_.endo_names = char(M_.endo_names, 'ID');
M_.endo_names_tex = char(M_.endo_names_tex, '{I_D}');
M_.endo_names_long = char(M_.endo_names_long, 'ID');
M_.endo_names = char(M_.endo_names, 'zeta');
M_.endo_names_tex = char(M_.endo_names_tex, '{\zeta}');
M_.endo_names_long = char(M_.endo_names_long, 'zeta');
M_.endo_names = char(M_.endo_names, 'RD');
M_.endo_names_tex = char(M_.endo_names_tex, '{\mathcal{R}_{D}}');
M_.endo_names_long = char(M_.endo_names_long, 'RD');
M_.endo_names = char(M_.endo_names, 'f_fcn');
M_.endo_names_tex = char(M_.endo_names_tex, '{\left.       f\left( \cdot \right)            \right|}');
M_.endo_names_long = char(M_.endo_names_long, 'f_fcn');
M_.endo_names = char(M_.endo_names, 'f_fcn_prime');
M_.endo_names_tex = char(M_.endo_names_tex, '{\left.       f^‎{\prime}\left( \cdot \right)   \right|}');
M_.endo_names_long = char(M_.endo_names_long, 'f_fcn_prime');
M_.endo_names = char(M_.endo_names, 'g_fcn');
M_.endo_names_tex = char(M_.endo_names_tex, '{\left.       g\left( \cdot \right)            \right|}');
M_.endo_names_long = char(M_.endo_names_long, 'g_fcn');
M_.endo_names = char(M_.endo_names, 'g_fcn_prime');
M_.endo_names_tex = char(M_.endo_names_tex, '{\left.       g^‎{\prime}\left( \cdot \right)   \right|}');
M_.endo_names_long = char(M_.endo_names_long, 'g_fcn_prime');
M_.param_names = 'beta';
M_.param_names_tex = '\beta';
M_.param_names_long = 'beta';
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, '\alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'epsilon');
M_.param_names_tex = char(M_.param_names_tex, '\epsilon');
M_.param_names_long = char(M_.param_names_long, 'epsilon');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, '\rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'varphi');
M_.param_names_tex = char(M_.param_names_tex, '\varphi');
M_.param_names_long = char(M_.param_names_long, 'varphi');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, '\delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'chi');
M_.param_names_tex = char(M_.param_names_tex, '\chi');
M_.param_names_long = char(M_.param_names_long, 'chi');
M_.param_names = char(M_.param_names, 'vartheta');
M_.param_names_tex = char(M_.param_names_tex, '\vartheta');
M_.param_names_long = char(M_.param_names_long, 'vartheta');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, '\gamma');
M_.param_names_long = char(M_.param_names_long, 'gamma');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, '\phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, '\eta');
M_.param_names_long = char(M_.param_names_long, 'eta');
M_.param_names = char(M_.param_names, 'LS');
M_.param_names_tex = char(M_.param_names_tex, 'LS');
M_.param_names_long = char(M_.param_names_long, 'LS');
M_.param_names = char(M_.param_names, 'lambda');
M_.param_names_tex = char(M_.param_names_tex, '\lambda');
M_.param_names_long = char(M_.param_names_long, 'lambda');
M_.param_names = char(M_.param_names, 'rhozeta');
M_.param_names_tex = char(M_.param_names_tex, '{\rho}_{\zeta}');
M_.param_names_long = char(M_.param_names_long, 'rhozeta');
M_.param_names = char(M_.param_names, 'M');
M_.param_names_tex = char(M_.param_names_tex, '\mathcal{M}');
M_.param_names_long = char(M_.param_names_long, 'M');
M_.param_names = char(M_.param_names, 'psi_N');
M_.param_names_tex = char(M_.param_names_tex, '\psi_N');
M_.param_names_long = char(M_.param_names_long, 'psi_N');
M_.param_names = char(M_.param_names, 'psi_I');
M_.param_names_tex = char(M_.param_names_tex, '\psi_I');
M_.param_names_long = char(M_.param_names_long, 'psi_I');
M_.param_names = char(M_.param_names, 'gg');
M_.param_names_tex = char(M_.param_names_tex, 'g^{BGP}');
M_.param_names_long = char(M_.param_names_long, 'gg');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 24;
M_.param_nbr = 18;
M_.orig_endo_nbr = 24;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('endogenous_growth_static');
erase_compiled_function('endogenous_growth_dynamic');
M_.lead_lag_incidence = [
 1 10 0;
 2 11 0;
 3 12 0;
 0 13 34;
 0 14 35;
 0 15 0;
 4 16 36;
 0 17 37;
 0 18 0;
 0 19 38;
 0 20 39;
 5 21 0;
 0 22 40;
 0 23 0;
 6 24 0;
 7 25 0;
 0 26 41;
 8 27 42;
 9 28 0;
 0 29 0;
 0 30 0;
 0 31 43;
 0 32 0;
 0 33 44;]';
M_.nstatic = 6;
M_.nfwrd   = 9;
M_.npred   = 7;
M_.nboth   = 2;
M_.nsfwrd   = 11;
M_.nspred   = 9;
M_.ndynamic   = 18;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(24, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(18, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 113;
M_.NNZDerivatives(2) = 245;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.96;
beta = M_.params( 1 );
M_.params( 2 ) = 0.33;
alpha = M_.params( 2 );
M_.params( 3 ) = 0.5;
epsilon = M_.params( 3 );
M_.params( 4 ) = 1;
rho = M_.params( 4 );
M_.params( 5 ) = 0.5;
varphi = M_.params( 5 );
M_.params( 5 ) = 5e-05;
varphi = M_.params( 5 );
M_.params( 6 ) = 0.10;
delta = M_.params( 6 );
M_.params( 7 ) = 0.2236;
chi = M_.params( 7 );
M_.params( 8 ) = 1+1/(1-M_.params(2));
vartheta = M_.params( 8 );
M_.params( 9 ) = 0.10;
gamma = M_.params( 9 );
M_.params( 9 ) = 0.12;
gamma = M_.params( 9 );
M_.params( 10 ) = 0.975;
phi = M_.params( 10 );
M_.params( 10 ) = 0.90;
phi = M_.params( 10 );
M_.params( 11 ) = 0.33;
eta = M_.params( 11 );
M_.params( 12 ) = 0.0645;
LS = M_.params( 12 );
M_.params( 13 ) = 0.1;
lambda = M_.params( 13 );
M_.params( 14 ) = 0.5;
rhozeta = M_.params( 14 );
M_.params( 15 ) = 1.315756236185665;
M = M_.params( 15 );
M_.params( 16 ) = 1;
psi_N = M_.params( 16 );
M_.params( 17 ) = 1;
psi_I = M_.params( 17 );
delete('ss_flexprice.mat');
disp('Steady State Solution:')
ss_solver_flex_price_model( M_.param_names, M_.params);
load ss_flexprice.mat
M_.params( 18 ) = steady_g;
gg = M_.params( 18 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = steady_g;
oo_.steady_state( 3 ) = steady_VD;
oo_.steady_state( 2 ) = steady_ZD;
oo_.steady_state( 19 ) = steady_zeta;
oo_.steady_state( 7 ) = steady_ND;
oo_.steady_state( 11 ) = steady_Lambda;
oo_.steady_state( 4 ) = steady_J;
oo_.steady_state( 5 ) = steady_H;
oo_.steady_state( 6 ) = steady_Pi;
oo_.steady_state( 8 ) = steady_YDW;
oo_.steady_state( 9 ) = steady_YD;
oo_.steady_state( 17 ) = steady_Q;
oo_.steady_state( 16 ) = steady_KD;
oo_.steady_state( 14 ) = steady_L;
oo_.steady_state( 18 ) = steady_ID;
oo_.steady_state( 10 ) = steady_CD;
oo_.steady_state( 15 ) = steady_GammaD;
oo_.steady_state( 13 ) = steady_muD;
oo_.steady_state( 12 ) = steady_UCD;
oo_.steady_state( 21 ) = 0;
oo_.steady_state( 23 ) = 0;
oo_.steady_state( 22 ) = 0;
oo_.steady_state( 24 ) = 0;
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.sigma_e_is_diagonal = 1;
options_.solve_algo = 1;
options_.steady.maxit = 50000;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
set_dynare_seed(092677);
save('endogenous_growth_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('endogenous_growth_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('endogenous_growth_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('endogenous_growth_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('endogenous_growth_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
