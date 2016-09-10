
pvarcoirfs = pvarcoirfs_clean;

% Put VAR IRFs in the same format (and delete final obs)
pvar_irf_rd     =  pvarcoirfs(strmatch('rd : rd', pvarcoirfs.id1), :).irf / 100;
pvar_irf_gdp    =  pvarcoirfs(strmatch('rd : gdp', pvarcoirfs.id1), :).irf / 100;
pvar_irf_sp     =  pvarcoirfs(strmatch('rd : sp', pvarcoirfs.id1), :).irf / 100;
pvar_irf_tfp    =  pvarcoirfs(strmatch('rd : tfp', pvarcoirfs.id1), :).irf / 100;

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
   plot( pvar_irf_rd, 'k--', 'linewidth', 2); 
   hold off;
   axis tight;
   jj = jj + 1;
   
   % Set plot axes (set min to zero)
   ymin = min(min( pvar_irf_rd ), 0);
   ymax = max(pvar_irf_rd );
   ylim([ymin ymax])

   subplot(2,2,jj);
   hold on;
   plot( pvar_irf_tfp, 'k--', 'linewidth', 2); 
   hold off;
   axis tight;
   jj = jj + 1;

   subplot(2,2,jj);
   hold on;
   plot( pvar_irf_gdp, 'k--', 'linewidth', 2); 
   hold off;
   axis tight;
   jj = jj + 1;

   subplot(2,2,jj);
   hold on;
   plot( pvar_irf_sp, 'k--', 'linewidth', 2); 
   hold off;
   axis tight;
  jj = jj + 1;




