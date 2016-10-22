%% This comes from resid.m
% Need to be a separate script because resid.m will not return residuals if
% no steady state found

% output : z (residuals)

% Compute the residuals
if options_.block && ~options_.bytecode
    z = zeros(M_.endo_nbr,1);
    for i = 1:length(M_.block_structure_stat.block)
        [r, g, yy, var_indx] = feval([M_.fname '_static'],...
                                     i,...
                                     oo_.steady_state,...
                                     [oo_.exo_steady_state; ...
                            oo_.exo_det_steady_state], M_.params);
        idx = M_.block_structure_stat.block(i).equation;
        z(idx) = r;
    end
elseif options_.bytecode
    [check, z] = bytecode('evaluate','static');
    mexErrCheck('bytecode', check);
else
    z = feval([M_.fname '_static'],...
              oo_.steady_state,...
              [oo_.exo_steady_state; ...
               oo_.exo_det_steady_state], M_.params);
end

