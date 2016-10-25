
%=========================================================================%
%%%%                       SIMULATED SERIES                            %%%%
%=========================================================================%

% Clear out the loglinear dr
load level0workspace oo_ options_

% Produce simulated series WITHOUT the loglinear command
% stoch_simul(order=1,periods=50000, irf=31, nograph, nodisplay, nocorr, nomoments, noprint, irf_shocks = (epsilon_n));

% Run stoch_simul funtion from dynare
options_.nocorr = 1;
options_.nograph = 1;
options_.nomoments = 1;
options_.order = 1;

options_.periods = 600;
options_.periods = 5000;

options_.irf = 11;
options_.loglinear = 0;                      %%% TURN OFF LOGLINEAR
options_.nodisplay = 0;
options_.noprint = 0;
options_.irf_shocks=[];
options_.irf_shocks = 'epsilon_n';

steady;

var_list_=[];
info = stoch_simul(var_list_);


% Collect the key series
g_simul = oo_.endo_simul(1,:);
RD_simul = oo_.endo_simul(24,:);
YD_simul = oo_.endo_simul(11,:);

% Compute TFP (using a cumulative product of g)
A_lead = cumprod(g_simul);
A_simul = ones(size(A_lead));
A_simul(2:end) = A_lead(1:end-1);

% Multiply this trend with R_D and Y_D
R_simul = RD_simul .* A_simul;
Y_simul = YD_simul .* A_simul;

%=========================================================================%
%%%%                     PLOT SIMULATED SERIES                         %%%%
%=========================================================================%

% % Plot
% plot_length = 300;
% figure();
% subplot(2,2,1);
% plot(R_simul(:,1:plot_length));
% title('R');
% 
% subplot(2,2,2);
% plot(A_simul(:,1:plot_length));
% title('A');
% 
% subplot(2,2,3);
% plot(Y_simul(:,1:plot_length));
% title('Y');
% 
% Output to a csv
SimulData = [log(A_simul'), log(R_simul'), log(Y_simul')];
csvwrite('Data_A_R_Y_logs.csv', SimulData);

% var_given_model(A_simul, R_simul, Y_simul)

% Sanity check: these are all equal to the steady state values
% mean(YD_simul);
% mean(g_simul);
% mean(RD_simul);