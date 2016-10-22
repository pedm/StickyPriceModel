function [ fd ] = first_difference( simul )
%Apply first difference. Requires Tx1 input, gives T-1x1 output

fd = simul(2:end) - simul(1:end-1);

end

