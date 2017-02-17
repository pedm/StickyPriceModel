%% Load PVAR IRFs

global pvarcoirfs_clean;

% Original panel var results (cum irfs of the var in first differences) Oct 2016
% load 'pvar_coirfs_full'                         % has 10 periods of data. pvar in first differences with gdp first
% load 'pvar_coirfs_full_30periods_tfp_first' 	% has 30 periods of data. pvar in first differences with tfp first
% pvarcoirfs_clean = pvarcoirfs;

% Panel VAR in levels Feb 2017
load 'PanelVAR_IRFs_Feb16_2017';             % has 10 periods of data. panel var in levels
pvarcoirfs_clean = PanelVAR_IRFs;            % it's not actually a cum irf

%% Plot PVAR IRFs
pvarcoirfs = pvarcoirfs_clean;

% drop all IRF steps beyond irf_length
% to modify this, go to estimation_init.m
try
    irf_length = options_.EST.irf_length;
catch
    options_.EST = [];
    options_.EST.irf_length = 11;
    irf_length = 11;
end

pvarcoirfs(pvarcoirfs.step >= irf_length, :) = [];

% Put VAR IRFs in the same format (and delete final obs)
pvar_irf_rd     =  pvarcoirfs(strmatch('rd : rd', pvarcoirfs.id1), :).irf / 100;
pvar_irf_gdp    =  pvarcoirfs(strmatch('rd : gdp', pvarcoirfs.id1), :).irf / 100;
pvar_irf_sp     =  pvarcoirfs(strmatch('rd : sp', pvarcoirfs.id1), :).irf / 100;
pvar_irf_tfp    =  pvarcoirfs(strmatch('rd : tfp', pvarcoirfs.id1), :).irf / 100;

pvar_irf_rd_se  = pvarcoirfs(strmatch('rd : rd', pvarcoirfs.id1), :).se /100;
pvar_irf_gdp_se = pvarcoirfs(strmatch('rd : gdp', pvarcoirfs.id1), :).se/100;
pvar_irf_sp_se  = pvarcoirfs(strmatch('rd : sp', pvarcoirfs.id1), :).se/100;
pvar_irf_tfp_se = pvarcoirfs(strmatch('rd : tfp', pvarcoirfs.id1), :).se/100;

% Plot the important variables
try
    h = findobj('name', 'important variables');
    figure(h);
catch
    figure('name', 'important variables');
end
jj = 1;
   subplot(2,2,jj);
   hold on;
   plot( pvar_irf_rd + pvar_irf_rd_se, 'k:', 'linewidth', 1.66); 
   plot( pvar_irf_rd - pvar_irf_rd_se, 'k:', 'linewidth', 1.66);
   plot( pvar_irf_rd, 'k--', 'linewidth', 2); hold on;

   hold off;
   axis tight;
   jj = jj + 1;
   
   % Set plot axes (set min to zero)
%    ymin = min(min( pvar_irf_rd ), 0);
%    ymax = max(pvar_irf_rd );
%    ylim([ymin ymax])

   subplot(2,2,jj);
   hold on;
   plot( pvar_irf_tfp + pvar_irf_tfp_se, 'k:', 'linewidth', 1.66);
   plot( pvar_irf_tfp - pvar_irf_tfp_se, 'k:', 'linewidth', 1.66);
   plot( pvar_irf_tfp, 'k--', 'linewidth', 2); 

   hold off;
   axis tight;
   jj = jj + 1;

   subplot(2,2,jj);
   hold on;
   plot( pvar_irf_gdp + pvar_irf_gdp_se, 'k:', 'linewidth', 1.66); 
   plot( pvar_irf_gdp - pvar_irf_gdp_se, 'k:', 'linewidth', 1.66);
   plot( pvar_irf_gdp, 'k--', 'linewidth', 2); 

   hold off;
   axis tight;
   jj = jj + 1;

   subplot(2,2,jj);
   hold on;  
   plot( pvar_irf_sp + pvar_irf_sp_se, 'k:', 'linewidth', 1.66); 
   plot( pvar_irf_sp - pvar_irf_sp_se, 'k:', 'linewidth', 1.66);
   plot( pvar_irf_sp, 'k--', 'linewidth', 2); 
   hold off;
   axis tight;
  jj = jj + 1;




