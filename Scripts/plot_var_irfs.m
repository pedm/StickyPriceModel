
pvarcoirfs = pvarcoirfs_clean;

% drop all IRF steps beyond irf_length
% to modify this, go to estimation_init.m
irf_length = options_.EST.irf_length;
pvarcoirfs(pvarcoirfs.step >= irf_length, :) = [];

% Put VAR IRFs in the same format (and delete final obs)
pvar_irf_rd     =  pvarcoirfs(strmatch('rd : rd', pvarcoirfs.id1), :).irf / 100;
pvar_irf_gdp    =  pvarcoirfs(strmatch('rd : gdp', pvarcoirfs.id1), :).irf / 100;
pvar_irf_sp     =  pvarcoirfs(strmatch('rd : sp', pvarcoirfs.id1), :).irf / 100;
pvar_irf_tfp    =  pvarcoirfs(strmatch('rd : tfp', pvarcoirfs.id1), :).irf / 100;

pvar_irf_rd_se  = pvarcoirfs(strmatch('rd : rd', pvarcoirfs.id1), :).se;
pvar_irf_gdp_se = pvarcoirfs(strmatch('rd : gdp', pvarcoirfs.id1), :).se;
pvar_irf_sp_se  = pvarcoirfs(strmatch('rd : sp', pvarcoirfs.id1), :).se;
pvar_irf_tfp_se = pvarcoirfs(strmatch('rd : tfp', pvarcoirfs.id1), :).se;

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




