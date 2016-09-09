
%===================================================================%
%%%%            COMPARE MODEL AND VAR IRFs                       %%%%
%===================================================================%

% disp('POST PROCESSING IRF DISTANCE:')

% Put VAR IRFs in the same format (and delete final obs)
% load('pvar_coirfs.mat')
% pvar_irf_rd = pvar_irf_rd(1:end-1)'/100;
% pvar_irf_gdp = pvar_irf_gdp(1:end-1)'/100;
% pvar_irf_sp = pvar_irf_sp(1:end-1)'/100;
% pvar_irf_tfp = pvar_irf_tfp(1:end-1)'/100;
% save 'pvar_coirfs_clean' pvar_irf_*

% load R&D IRFs from VAR
load 'pvar_coirfs_clean';

% Select the model irfs
ii = find(ismember(IRnames_dynare, 'R'));
mod_irf_rd = IR_dynare(ii,:);

ii = find(ismember(IRnames_dynare, 'A'));
mod_irf_tfp = IR_dynare(ii,:);

ii = find(ismember(IRnames_dynare, 'Y'));
mod_irf_gdp = IR_dynare(ii,:);

ii = find(ismember(IRnames_dynare, 'S'));
mod_irf_sp = IR_dynare(ii,:);

% Calculate the distance between irfs
irf_distance = sum(abs(mod_irf_rd - pvar_irf_rd)) + sum(abs(mod_irf_tfp - pvar_irf_tfp)) + sum(abs(mod_irf_gdp - pvar_irf_gdp)) + sum(abs(mod_irf_sp - pvar_irf_sp))

% TODO: plot a comparison between model and var irfs

% TODO: perhaps I should use quadratic distance
