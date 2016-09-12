
%===================================================================%
%%%%            COMPARE MODEL AND VAR IRFs                       %%%%
%===================================================================%

% create a local copy of VAR IRFs (pvarcoirfs_clean is a global, so I only
% have to load it once)
pvarcoirfs = pvarcoirfs_clean;

% Select the model irfs
ii = find(ismember(IRnames_dynare, 'R'));
mod_irf_rd = IR_dynare(ii,:);

ii = find(ismember(IRnames_dynare, 'A'));
mod_irf_tfp = IR_dynare(ii,:);

ii = find(ismember(IRnames_dynare, 'Y'));
mod_irf_gdp = IR_dynare(ii,:);

ii = find(ismember(IRnames_dynare, 'S'));
mod_irf_sp = IR_dynare(ii,:);

%% Add model IRFs to table

pvarcoirfs(strmatch('rd : sp', pvarcoirfs.id1), :).model    = mod_irf_sp';
pvarcoirfs(strmatch('rd : rd', pvarcoirfs.id1), :).model    = mod_irf_rd';
pvarcoirfs(strmatch('rd : tfp', pvarcoirfs.id1), :).model   = mod_irf_tfp';
pvarcoirfs(strmatch('rd : gdp', pvarcoirfs.id1), :).model   = mod_irf_gdp';

% Scale by 100
pvarcoirfs.model = pvarcoirfs.model*100;

% drop table rows with 0 standard errors
pvarcoirfs(pvarcoirfs.se == 0, :) = [];

% Calculate the distance between irfs
% Use same formula as CEE
DDD = pvarcoirfs.irf - pvarcoirfs.model;
VVV = diag(pvarcoirfs.se.*pvarcoirfs.se); % create diagonal matrix of the variances
irf_distance_sub = DDD'* inv(VVV) * DDD;

%% Add growth rate
ss = exp(oo_.steady_state);
growth_rate = ss(1);

irf_distance = irf_distance_sub + abs(growth_rate - 1.0118)*10e+08;

% TODO: add more weight to the growth rate
