function [ h_simul ] = helmert( simul )
% Apply Helmert's Transformation from Arellano and Bover 1995
% Input must be Tx1. Output will be (T-1)x1

% http://www.cemfi.es/~arellano/arellano-bover-1995.pdf

T = length(simul);
c_list = 1:T-1;
c = (T - c_list)./(T - c_list + 1);

% diag T-1 x T-1 matrix
diag_ones = diag(ones(1,T-1));

% add extra column
diag_ones(:,T) = 0;

% create upper partial matrix
A_part = triu((-(T-c_list').^-1)*ones(1,T));

% replace diag with ones
A_part(diag_ones == 1) = 1;

% check that all rows sum to one
% sum(A_part,2)

% Combine
A = diag(sqrt(c)) * A_part;

% Apply the transformation
h_simul = A*simul;

end

