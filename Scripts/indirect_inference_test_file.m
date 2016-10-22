%% Transform Series

log_R = log( R_simul(1:20)' );
log_R(3) = NaN

% replace any cases where R_simul < 0

fd_R = first_difference( log_R );
h_R = helmert(fd_R);


% Load
% [xlsdata, xlstext] = xlsread('SW2001_Data.xls','Sheet1');
% Define data
X = [log(A_simul'), log(R_simul'), log(Y_simul')];
% Define label for plots
% dates = xlstext(2:end,1);
% vnames = xlstext(1,2:end);
% Define number of variables and of observations
[nobs, nvar] = size(X);


% TODO: still working on this


%% VAR ESTIMATION
% =======================================================================
% Set the case for the VARout (0, 1, or 2)
det = 2;
% Set number of nlags
nlags = 4;
% Estimate 
[VAR, VARopt] = VARmodel(X,nlags,det);
% Print at screen and create table
VARopt.vnames = vnames;
[beta, tstat, TABLE] = VARprint(VAR,VARopt);


%% IMPULSE RESPONSE
% =======================================================================
% Set options some options for IRF calculation
VARopt.nsteps = 60;
VARopt.ident = 'oir';
VARopt.quality = 0;
% Compute IRF
[IRF, VAR] = VARir(VAR,VARopt);
% Compute error bands
[IRFINF,IRFSUP,IRFMED] = VARirband(VAR,VARopt);
% Plot
VARirplot(IRFMED,VARopt,IRFINF,IRFSUP);
