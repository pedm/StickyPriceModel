clc; clear all; close all;


% set model parameters ----------------------------------------------------

% conventional

param.beta   = 0.96; 
param.alpha  = 0.33;
param.epsilon= 1/2; % Inverse Frisch labor supply elasticity
param.rho    =   1; % Inverse of intertemporal elasticity of substitution
param.varphi = 0.5 / 10000; % elasticity of Q to I/K ratio 
param.delta  = 0.10;
param.chi    = 0.05^param.epsilon; % disutility of labor supply
param.chi    = 0.2236; % dynare param
param.vartheta = 1 + 1/(1-param.alpha) ;

% growth

param.gamma  = 0.12; % elasticity of labor disutility to technology
% param.gamma  = 0.10;  %dynare param
param.phi    = 0.90;
param.eta    = 0.33;
param.LS     = 0.0645; %0.07; 
param.lambda = .25; % adoption probability
param.lambda = 0.1; %dynare param
param.M      = 4.167 / (4.167 - 1);
% shocks

param.rhozeta = 0.50; % persistence of exogenous "innovation productivity" shock

% -------------------------------------------------------------------------


% STEADY STATE
ss = steady(param);
disp('Steady State Solution:')
ss

