function [ value ] = param_value( param_name )

global M_
value = M_.params(strcmp(param_name, cellstr(M_.param_names)));
disp(M_.param_names(strcmp(param_name, cellstr(M_.param_names)), :));

end

