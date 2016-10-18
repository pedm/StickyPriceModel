
% close all;
% disp('POST PROCESSING IRFs')

%===================================================================%
%%%%            PRODUCE CUMMULATIVE IRFS                         %%%%
%===================================================================%

% figure('name', 'cummulative irfs')

% In linear form, detrended variables are defined as:
% Y_D = Y / A
% And growth rate g equals:
% g = A(+1) / A
% A(+1) = g * A

% With loglinear option specified, everything becomes: 
% logY = logY_D + logA
% logA = log(g(-1)) + logA(-1)

% Compute productivity (lagged cumsum of productivity growth log(g) )
A_lead = cumsum(oo_.irfs.g_epsilon_chi);
A = zeros(size(A_lead));
A(2:end) = A_lead(1:end-1);
oo_.irfs.A_epsilon_chi = A;

% Modify the detrended variables
oo_.irfs.Z_epsilon_chi = oo_.irfs.ZD_epsilon_chi + A;
oo_.irfs.V_epsilon_chi = oo_.irfs.VD_epsilon_chi + A;
oo_.irfs.N_epsilon_chi = oo_.irfs.ND_epsilon_chi + A;
oo_.irfs.Y_epsilon_chi = oo_.irfs.YD_epsilon_chi + A;
oo_.irfs.C_epsilon_chi = oo_.irfs.CD_epsilon_chi + A;
oo_.irfs.K_epsilon_chi = oo_.irfs.KD_epsilon_chi + A;
oo_.irfs.I_epsilon_chi = oo_.irfs.ID_epsilon_chi + A;
oo_.irfs.S_epsilon_chi = oo_.irfs.SD_epsilon_chi + A;
oo_.irfs.X_epsilon_chi = oo_.irfs.XD_epsilon_chi + A;
oo_.irfs.R_epsilon_chi = oo_.irfs.RD_epsilon_chi + A;
oo_.irfs.UC_epsilon_chi = oo_.irfs.UCD_epsilon_chi + A;
oo_.irfs.mu_epsilon_chi = oo_.irfs.muD_epsilon_chi + A;
oo_.irfs.Gamma_epsilon_chi = oo_.irfs.GammaD_epsilon_chi + A;

% Prep Data so I can plot it in same format as Albert's main.m

% Save names
IRnames_dynare = fieldnames(oo_.irfs);

% Save data
IR_dynare = NaN(length(IRnames_dynare), options_.EST.irf_length);

for ii = 1:length(IRnames_dynare)
    full_irf = getfield(oo_.irfs, char(IRnames_dynare(ii)));
    IR_dynare(ii, :) = full_irf(1:options_.EST.irf_length);
end
IRnames_dynare = strrep(IRnames_dynare, '_epsilon_chi', '');
IRnames_dynare = strrep(IRnames_dynare, '_', '');
save IR_dynare.mat IR_dynare IRnames_dynare


