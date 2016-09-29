LB.eta = 0;  
LB.gamma = 0.00001;
LB.phi = 0.5;
LB.lambda_bar = 0.005;
LB.psi_N = 0;
LB.rhozeta = 0.0001; 
LB.rhozeta2 = 0.0001;  
LB.sigmazeta = .5;
LB.zetabar = .1;

UB.eta = 0.9999;  
UB.gamma = 0.9999; 
UB.phi = 0.99; 
UB.lambda_bar = 10; 
UB.psi_N = 100; 
UB.rhozeta = 1.5; 
UB.rhozeta2 = 0.99; 
UB.sigmazeta = 10;  
UB.zetabar = 10;

FF = fields(LB)
 
DIST.(char(FF(1))) = UB.(char(FF(1))) - LB.(char(FF(1)));
DIST.(char(FF(2))) = UB.(char(FF(2))) - LB.(char(FF(2)));
DIST.(char(FF(3))) = UB.(char(FF(3))) - LB.(char(FF(3)));
DIST.(char(FF(4))) = UB.(char(FF(4))) - LB.(char(FF(4)));
DIST.(char(FF(5))) = UB.(char(FF(5))) - LB.(char(FF(5)));
DIST.(char(FF(6))) = UB.(char(FF(6))) - LB.(char(FF(6)));
DIST.(char(FF(7))) = UB.(char(FF(7))) - LB.(char(FF(7)));
DIST.(char(FF(8))) = UB.(char(FF(8))) - LB.(char(FF(8)));
DIST.(char(FF(9))) = UB.(char(FF(9))) - LB.(char(FF(9)));

v={
	LB.(char(FF(1))) + 0.2 * DIST.(char(FF(1))) : 0.6 * DIST.(char(FF(1))) : UB.(char(FF(1))) - 0.2 * DIST.(char(FF(1)))
	LB.(char(FF(2))) + 0.2 * DIST.(char(FF(2))) : 0.6 * DIST.(char(FF(2))) : UB.(char(FF(2))) - 0.2 * DIST.(char(FF(2)))
	LB.(char(FF(3))) + 0.2 * DIST.(char(FF(3))) : 0.6 * DIST.(char(FF(3))) : UB.(char(FF(3))) - 0.2 * DIST.(char(FF(3)))
	LB.(char(FF(4))) + 0.2 * DIST.(char(FF(4))) : 0.6 * DIST.(char(FF(4))) : UB.(char(FF(4))) - 0.2 * DIST.(char(FF(4)))
	LB.(char(FF(5))) + 0.2 * DIST.(char(FF(5))) : 0.6 * DIST.(char(FF(5))) : UB.(char(FF(5))) - 0.2 * DIST.(char(FF(5)))
	LB.(char(FF(6))) + 0.2 * DIST.(char(FF(6))) : 0.6 * DIST.(char(FF(6))) : UB.(char(FF(6))) - 0.2 * DIST.(char(FF(6)))
	LB.(char(FF(7))) + 0.2 * DIST.(char(FF(7))) : 0.6 * DIST.(char(FF(7))) : UB.(char(FF(7))) - 0.2 * DIST.(char(FF(7)))
	LB.(char(FF(8))) + 0.2 * DIST.(char(FF(8))) : 0.6 * DIST.(char(FF(8))) : UB.(char(FF(8))) - 0.2 * DIST.(char(FF(8)))
	LB.(char(FF(9))) + 0.2 * DIST.(char(FF(9))) : 0.6 * DIST.(char(FF(9))) : UB.(char(FF(9))) - 0.2 * DIST.(char(FF(9)))
};



 
% ii = 3
% percent = 0.2
% distance = UB.(char(FF(ii))) - LB.(char(FF(ii)));
% pp = LB.(char(FF(ii))) + percent * distance


%% Permute
% the data (slightly different)
%      v={
%           1:2
%           1:2
%           1:2
%           1:2
%      };
% the engine
     n_count=numel(v);
     x_cell=cell(n_count,1);
     [x_cell{1:n_count,1}]=ndgrid(v{end:-1:1});
     permutations=reshape(cat(n_count+1,x_cell{:}),[],n_count);
     permutations=permutations(:,end:-1:1);
% the result
     % disp(permutations);

     
%% Guess


for kk = 1:length(permutations);
    guess = permutations(kk,:);

    set_param_value('eta', guess(1) );
    set_param_value('gamma', guess(2) );
    set_param_value('phi', guess(3) );
    set_param_value('lambda_bar', guess(4) );
    set_param_value('psi_N', guess(5) );
    set_param_value('rhozeta', guess(6) );
    set_param_value('rhozeta2', guess(7) );
    set_param_value('sigmazeta', guess(8) );
    set_param_value('zetabar', guess(9) );

    try
        var_list_=[];
        info = stoch_simul(var_list_);

        disp('Success!')
        guess
        kk
        return
    catch
        disp('No ss found')
        kk
    end
    

end

