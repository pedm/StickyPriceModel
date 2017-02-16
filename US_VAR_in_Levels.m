%% Setup 
clear;clc;close all;
addpath('VAR')
addpath('VAR\VAR')
addpath('VAR\Utils')
addpath('VAR\Stats')
addpath('VAR\Auxiliary')
addpath('VAR\Figure')
addpath('VAR\ExportFig')

%% US Data
Data = readtable('US_logged_series.csv');
X = [Data.gdp, Data.tfp, Data.rd];
vnames = {'GDP', 'TFP', 'R&D'};

% Define number of variables and of observations
[nobs, nvar] = size(X);

%% VAR ESTIMATION
% =======================================================================
% Set the case for the VARout (0, 1, or 2)
% const: 0 no constant; 1 constant; 2 constant and trend; 3 constant, 
%       trend, and trend^2 [dflt = 0]
const = 1;
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

%% Compute error bands
% VARopt.pctg = 95;
[IRFINF,IRFSUP,IRFMED] = VARirband(VAR,VARopt);

%% Plot
VARirplot(IRF,VARopt, IRFINF, IRFSUP);
