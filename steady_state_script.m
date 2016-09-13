clear; clc;
load('endogenous_growth_results.mat');
M_temp = M_;
clear M_;

%% Set Params
param.beta    = 0.96; 
param.alpha   = 0.33;
param.epsilon = 1/2;                   % Inverse Frisch labor supply elasticity
param.rho     =   1;                   % Inverse of intertemporal elasticity of substitution
param.delta   = 0.10;
param.chi     = 1.5652;                % Disutility of labor supply
param.vartheta = 1 + 1/(1-param.alpha);
param.eta    = 0.999;       % 0.375;                   % Curvature of innovations production in R&D expenditure (original = 0.33)
param.gamma  = 0.0396;       % 0.35;                    % weight of current consumption on JR term; indexes strength of wealth effects (0->no wealth effect (GHH), 1-> KPR prefs)
param.phi    = 0.9862;       % 0.875;                   % Survival rate of technologies
param.lambda = 0.4466;       % 0.075;                   % Adoption probability
param.rhozeta    = 0.5487; % 0.5; 
param.rhozeta2   = 0.0004; % 0.1;                 % Note: there's a minus in front of this (also, in estimation, must be greater than 0)
param.sigmazeta  = 1.0613; % 3.5;
param.zetabar    = 0.4777;
param.M       = 4.167 / (4.167 - 1);         % Markup. In the flex price model, markup is exogenous and given by M = ω/(ω − 1). I took this numbers from Gertler-Karadi “a model of unconventional monetary policy�?, who take them from estimates by Primiceri et al
param.psi_N   = 21.8893;                           % Adjustment cost to N
param.psi_I   = 1;                           % Adjustment cost to I

%% Define M_ object
M_ = [];
M_.param_names = fieldnames(param);
M_.params = struct2array(param)';
M_.param_nbr = length(M_.params);
M_.endo_names = M_temp.endo_names;

%% Compute SS
steady = endogenous_growth_steadystate();

ss = [cellstr(M_.endo_names), num2cell(steady)]
