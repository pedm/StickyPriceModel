%% Transform Series

tic
R_sim = transform_series( R_simul' );
A_sim = transform_series( A_simul' );
Y_sim = transform_series( Y_simul' );

% Define data
X = [R_sim, Y_sim, A_sim];
vnames = {'R&D', 'GDP', 'TFP'};

% Define number of variables and of observations
[nobs, nvar] = size(X);

%% VAR ESTIMATION
% =======================================================================
% Set the case for the VARout (0, 1, or 2)
% const: 0 no constant; 1 constant; 2 constant and trend; 3 constant, 
%       trend, and trend^2 [dflt = 0]
const = 0;
% Set number of nlags
nlags = 1;

% Estimate 
% NOTE: I added list wise deletion to VARmodel.m (drop any rows with missing data)
[VAR, VARopt] = VARmodel(X,nlags,const);
% Print at screen and create table
VARopt.vnames = vnames;
[beta, tstat, TABLE] = VARprint(VAR,VARopt);


%% IMPULSE RESPONSE
% =======================================================================
% Set options some options for IRF calculation
VARopt.nsteps = 10;
VARopt.ident = 'oir';
VARopt.quality = 0;
% Compute IRF
[IRF, VAR] = VARir(VAR,VARopt);

toc


% Compute error bands
% [IRFINF,IRFSUP,IRFMED] = VARirband(VAR,VARopt);
% Plot
VARirplot(IRF,VARopt);
