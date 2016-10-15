%% if loglinear option is selected
% g_simul = oo_.endo_simul(1,:)
% A_lead = cumsum(g_simul)
% A = zeros(size(A_lead));
% A_simul(2:end) = A_lead(1:end-1);
% A_simul = zeros(size(A_lead));
% A_simul(2:end) = A_lead(1:end-1);
% A_simul
% plot(A_simul)
% close all
% plot(A_simul)
% M_.endo_names
% M_.endo_names(11,:)
% M_.endo_names(24,:)
% M_.endo_names(22,:)
% RD_simul = oo_.endo_simul(24,:)
% size(RD_simul)
% size(A_simul)
% R_simul = RD_simul + A_simul
% plot(R_simul)
% hold on
% plot(RD_simul)
% hold off
% plot(exp(R_simul))
% plot(RD_simul)
% plot(R_simul)
% oo_.steady_state
% oo_.steady_state(24)
% exp(oo_.steady_state(24))
% plot(exp(A_simul))

%% if not
close all

g_simul = oo_.endo_simul(1,:);
A_lead = cumprod(g_simul);
A_simul = ones(size(A_lead));
A_simul(2:end) = A_lead(1:end-1)

RD_simul = oo_.endo_simul(24,:)
R_simul = RD_simul .* A_simul


YD_simul = oo_.endo_simul(11,:)
Y_simul = YD_simul .* A_simul

% Yup, looks good. it's revolving around steady state
mean(YD_simul)
oo_.steady_state(11)

plot(Y_simul(:,1:100))

Data = [A_simul', R_simul', Y_simul']
csvwrite('Data_A_R_Y.csv', Data)

