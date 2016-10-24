function [ h_simul ] = helmert( simul )
% Apply Helmert's Transformation from Arellano and Bover 1995
% Input must be Tx1. Output will be (T-1)x1
% http://www.cemfi.es/~arellano/arellano-bover-1995.pdf

T = length(simul);
oo = ones(size(simul));
oo(isnan(simul)) = 0;

t_list = 1:T-1;
c = (T - t_list)./(T - t_list + 1);

%% Cumsum Method (faster)
n = cumsum(oo, 'reverse') - 1;
n(n<=0) = NaN;
simul_nonan = simul;
simul_nonan(isnan(simul)) = 0;
csum = cumsum(simul_nonan, 'reverse');
% h_simul = sqrt(c)' .* ( simul(1:end-1) - ((1./(T - t_list')).*csum(2:end)) );
m = (csum-simul) ./ n;
w = sqrt(n ./ (n+1));

h_simul_ = w.*(simul - m);
h_simul = h_simul_(1:end-1);

%% Matrix Method (my original attempt. too slow)
% % diag T-1 x T-1 matrix
% diag_ones = diag(ones(1,T-1));
% 
% % add extra column
% diag_ones(:,T) = 0;
% 
% % create upper partial matrix
% A_part = triu((-(T-t_list').^-1)*ones(1,T));
% 
% % replace diag with ones
% A_part(diag_ones == 1) = 1;
% 
% % check that all rows sum to one
% % sum(A_part,2)
% 
% % Combine
% A = diag(sqrt(c)) * A_part;
% 
% % Apply the transformation
% h_simul2 = A*simul

end

