% This example shows how to compute IRFs, HDs, and FEVDs in a VAR with 
% data for inflation, unemployment, and interest rates.  Identification 
% is achieved by imposing short-run restrictions, computed with a Cholesky 
% decomposition of the reduced-form residuals' covariance matrix.  

% The VAR Toolbox 2.0 is required to run this code. To get the 
% latest version of the toolboxes visit: 
% 
%       https://sites.google.com/site/ambropo/MatlabCodes
% 
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com


%% 1. PRELIMINARIES
% =======================================================================
clear all; clear session; close all; clc
warning off all

% Load
[xlsdata, xlstext] = xlsread('SW2001_Data.xls','Sheet1');
% Define data
X = xlsdata;
% Define label for plots
dates = xlstext(2:end,1);
vnames = xlstext(1,2:end);
% Define number of variables and of observations
[nobs, nvar] = size(X);


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


%% HISTORICAL DECOMPOSITION
% =======================================================================
% Compute HD
HD = VARhd(VAR);
% Plot HD
VARhdplot(HD,VARopt);


%% FORECAST ERROR VARIANCE DECOMPOSITION
% =======================================================================
% Compute FEVD
[FEVD, VAR] = VARfevd(VAR,VARopt);
% Compute error bands
[FEVDINF,FEVDSUP,FEVDMED] = VARfevdband(VAR,VARopt);
% Plot
VARfevdplot(FEVDMED,VARopt,FEVDINF,FEVDSUP);

