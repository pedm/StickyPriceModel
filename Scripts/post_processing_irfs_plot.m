
%===================================================================%
%%%%            PLOT IRFS                                        %%%%
%===================================================================%

% Plot in same format as Albert
plot_cutoff = 36;
try
    h = findobj('name', 'original vars');
    figure(h);
catch
    figure('name', 'original vars');
    suptitle('Unmodified Dynare IRFs (Detrended)');
end
for ii = 1:plot_cutoff;
    subplot(6,6,ii);
    hold on;
    plot( IR_dynare(ii,:), 'linewidth', 2); 
    hold off;
    title( IRnames_dynare(ii) );
    axis tight;
end

% % Plot in same format as Albert
% figure('name', 'all variables dynare');
% for ii = 1:size(IR_dynare,1);
%     subplot(7,6,ii);
%     plot( IR_dynare(ii,:), 'linewidth', 2); 
%     title(IRnames_dynare(ii) );
%     axis tight;
% end
% suptitle('Dynare IRFs (Unmodified)');


% Undetrended vars
try
    h = findobj('name', 'un detrended vars');
    figure(h);
catch
    figure('name', 'un detrended vars');
    suptitle('Modified IRFs (With Trend Added)');
end
for ii = plot_cutoff+1:size(IR_dynare,1);
    subplot(4,4,ii-plot_cutoff);
    hold on;
    plot( IR_dynare(ii,:), 'linewidth', 2); 
    hold off;
    title(IRnames_dynare(ii) );
    axis tight;
end


% Plot the important variables
try
    h = findobj('name', 'important variables');
    figure(h);
catch
    figure('name', 'important variables');
end
jj = 1;
ii = find(ismember(IRnames_dynare, 'R'));
   subplot(2,2,jj);
   hold on;
   plot( IR_dynare(ii,:), 'linewidth', 2); 
   hold off;
   title(IRnames_dynare(ii) );
   axis tight;
   jj = jj + 1;
   
   % Set plot axes (set min to zero)
%    ymin = min(min(IR_dynare(ii,:) ), 0);
%    ymax = max(IR_dynare(ii,:) );
%    ylim([ymin ymax])

ii = find(ismember(IRnames_dynare, 'A'));
   subplot(2,2,jj);
   hold on;
   plot( IR_dynare(ii,:), 'linewidth', 2); 
   hold off;
   title(IRnames_dynare(ii) );
   axis tight;
   jj = jj + 1;

ii = find(ismember(IRnames_dynare, 'Y'));
   subplot(2,2,jj);
   hold on;
   plot( IR_dynare(ii,:), 'linewidth', 2); 
   hold off;
   title(IRnames_dynare(ii) );
   axis tight;
   jj = jj + 1;

ii = find(ismember(IRnames_dynare, 'S'));
   subplot(2,2,jj);
   hold on;
   plot( IR_dynare(ii,:), 'linewidth', 2); 
   hold off;
   title(IRnames_dynare(ii) );
   axis tight;
  jj = jj + 1;




