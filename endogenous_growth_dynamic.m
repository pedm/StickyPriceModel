function [residual, g1, g2, g3] = endogenous_growth_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           columns: equations in order of declaration
%                                                           rows: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(26, 1);
T59 = 1/params(7)*1/params(15);
T67 = (y(11)/y(16))^(1-params(10));
T85 = (y(10)*y(39)/y(16))^2;
T94 = (y(7)/y(1))^params(2);
T97 = y(23)^(1-params(2));
T117 = y(1)^(-params(4));
T124 = params(6)/(1+params(3));
T126 = y(23)^(1+params(3));
T128 = y(19)-y(24)*T124*T126;
T129 = T128^(-params(4));
T138 = (y(6)/(y(1)*y(19)))^(1-params(8));
T145 = y(10)^(-params(4))*y(43);
T149 = (y(10)*y(41)/y(24))^params(8);
T160 = y(24)*params(6)*y(23)^params(3)*1/y(21);
T165 = (1-params(2))*1/params(15)*(params(7)-1)/params(7);
T171 = (y(6)/y(1))^(1-params(8));
T191 = (y(10)*y(45)/y(27))^2;
T238 = params(16)/2;
T243 = exp(T238*(y(1)*y(16)/y(4)-params(18))^2);
T246 = exp(params(16)*(y(1)*y(16)/y(4)-params(18)));
T249 = params(17)/2;
T253 = exp(T249*(y(1)*y(27)/y(8)-params(18))^2);
T256 = exp(params(17)*(y(1)*y(27)/y(8)-params(18)));
T264 = getPowerDeriv(y(7)/y(1),params(2),1);
T275 = getPowerDeriv(y(6)/(y(1)*y(19)),1-params(8),1);
T281 = getPowerDeriv(y(6)/y(1),1-params(8),1);
T290 = 2*(y(1)*y(16)/y(4)-params(18));
T298 = 2*(y(1)*y(27)/y(8)-params(18));
T314 = getPowerDeriv(y(10)*y(41)/y(24),params(8),1);
T342 = getPowerDeriv(y(11)/y(16),1-params(10),1);
T360 = (-(y(1)*y(16)))/(y(4)*y(4));
T405 = getPowerDeriv(T128,(-params(4)),1);
T460 = T405*(-(y(24)*T124*getPowerDeriv(y(23),1+params(3),1)));
T523 = (-(y(1)*y(27)))/(y(8)*y(8));
lhs =y(10);
rhs =params(9)*(1+params(11)*(y(11)+y(12)-1));
residual(1)= lhs-rhs;
lhs =y(11)*y(1);
rhs =params(9)*(y(2)+(1-params(11))*y(3));
residual(2)= lhs-rhs;
lhs =y(12);
rhs =params(14)*y(28)*y(11)^(1-params(10))*y(16)^params(10);
residual(3)= lhs-rhs;
lhs =y(13);
rhs =params(11)*y(14)+params(9)*(1-params(11))*y(42)*y(37);
residual(4)= lhs-rhs;
lhs =y(14);
rhs =y(15)+params(9)*y(42)*y(38);
residual(5)= lhs-rhs;
lhs =y(15);
rhs =T59*y(17);
residual(6)= lhs-rhs;
lhs =y(28)*y(13)*params(14)*params(10)*T67;
rhs =1+log(y(33))*y(1)*y(16)/y(4)+log(y(32))-y(42)*log(y(47))*T85;
residual(7)= lhs-rhs;
lhs =y(18);
rhs =y(17);
residual(8)= lhs-rhs;
lhs =y(17);
rhs =T94*T97;
residual(9)= lhs-rhs;
lhs =y(18);
rhs =y(16)+y(19)+(1+log(y(34)))*y(27);
residual(10)= lhs-rhs;
lhs =y(20);
rhs =params(1)*y(21)/y(5)*T117;
residual(11)= lhs-rhs;
lhs =y(21);
rhs =T129+(-y(22))*params(8)*T138;
residual(12)= lhs-rhs;
lhs =y(22);
rhs =params(1)*(1-params(8))*T145*T149+T126*T124*T129;
residual(13)= lhs-rhs;
lhs =T129*T160;
rhs =T165*y(18)/y(23);
residual(14)= lhs-rhs;
lhs =y(24);
rhs =y(19)^params(8)*T171;
residual(15)= lhs-rhs;
lhs =y(25);
rhs =y(27)+y(7)/y(1)*(1-params(5));
residual(16)= lhs-rhs;
lhs =y(26);
rhs =1+log(y(34))+y(1)*y(27)/y(8)*log(y(35))-y(42)*T191*log(y(48));
residual(17)= lhs-rhs;
lhs =y(26);
rhs =y(42)*(params(2)*y(10)*(params(7)-1)*y(40)/(params(7)*params(15)*y(25))+(1-params(5))*y(44));
residual(18)= lhs-rhs;
lhs =log(y(28));
rhs =params(12)*log(y(9))+params(13)*x(it_, 1);
residual(19)= lhs-rhs;
lhs =y(29);
rhs =y(14)+y(25)*y(26)+(y(11)+y(12)-1)*y(13)+y(30);
residual(20)= lhs-rhs;
lhs =y(30);
rhs =y(10)*y(42)*(y(37)*y(36)+y(46));
residual(21)= lhs-rhs;
lhs =y(31);
rhs =y(16);
residual(22)= lhs-rhs;
lhs =y(32);
rhs =T243;
residual(23)= lhs-rhs;
lhs =y(33);
rhs =T246;
residual(24)= lhs-rhs;
lhs =y(34);
rhs =T253;
residual(25)= lhs-rhs;
lhs =y(35);
rhs =T256;
residual(26)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(26, 49);

  %
  % Jacobian matrix
  %

  g1(1,10)=1;
  g1(1,11)=(-(params(9)*params(11)));
  g1(1,12)=(-(params(9)*params(11)));
  g1(2,1)=y(11);
  g1(2,2)=(-params(9));
  g1(2,11)=y(1);
  g1(2,3)=(-(params(9)*(1-params(11))));
  g1(3,11)=(-(y(16)^params(10)*params(14)*y(28)*getPowerDeriv(y(11),1-params(10),1)));
  g1(3,12)=1;
  g1(3,16)=(-(params(14)*y(28)*y(11)^(1-params(10))*getPowerDeriv(y(16),params(10),1)));
  g1(3,28)=(-(y(16)^params(10)*params(14)*y(11)^(1-params(10))));
  g1(4,13)=1;
  g1(4,37)=(-(params(9)*(1-params(11))*y(42)));
  g1(4,14)=(-params(11));
  g1(4,42)=(-(params(9)*(1-params(11))*y(37)));
  g1(5,14)=1;
  g1(5,38)=(-(params(9)*y(42)));
  g1(5,15)=(-1);
  g1(5,42)=(-(params(9)*y(38)));
  g1(6,15)=1;
  g1(6,17)=(-T59);
  g1(7,1)=(-(log(y(33))*y(16)/y(4)));
  g1(7,10)=y(42)*log(y(47))*y(39)/y(16)*2*y(10)*y(39)/y(16);
  g1(7,11)=y(28)*y(13)*params(14)*params(10)*1/y(16)*T342;
  g1(7,13)=T67*y(28)*params(14)*params(10);
  g1(7,4)=(-(log(y(33))*T360));
  g1(7,16)=y(28)*y(13)*params(14)*params(10)*T342*(-y(11))/(y(16)*y(16))-(log(y(33))*y(1)/y(4)-y(42)*log(y(47))*2*y(10)*y(39)/y(16)*(-(y(10)*y(39)))/(y(16)*y(16)));
  g1(7,39)=y(42)*log(y(47))*2*y(10)*y(39)/y(16)*y(10)/y(16);
  g1(7,42)=log(y(47))*T85;
  g1(7,28)=y(13)*params(14)*params(10)*T67;
  g1(7,32)=(-(1/y(32)));
  g1(7,33)=(-(y(1)*y(16)/y(4)*1/y(33)));
  g1(7,47)=T85*y(42)*1/y(47);
  g1(8,17)=(-1);
  g1(8,18)=1;
  g1(9,1)=(-(T97*(-y(7))/(y(1)*y(1))*T264));
  g1(9,17)=1;
  g1(9,23)=(-(T94*getPowerDeriv(y(23),1-params(2),1)));
  g1(9,7)=(-(T97*T264*1/y(1)));
  g1(10,16)=(-1);
  g1(10,18)=1;
  g1(10,19)=(-1);
  g1(10,27)=(-(1+log(y(34))));
  g1(10,34)=(-(y(27)*1/y(34)));
  g1(11,1)=(-(params(1)*y(21)/y(5)*getPowerDeriv(y(1),(-params(4)),1)));
  g1(11,20)=1;
  g1(11,5)=(-(T117*(-(params(1)*y(21)))/(y(5)*y(5))));
  g1(11,21)=(-(T117*params(1)/y(5)));
  g1(12,1)=(-((-y(22))*params(8)*(-(y(19)*y(6)))/(y(1)*y(19)*y(1)*y(19))*T275));
  g1(12,19)=(-(T405+(-y(22))*params(8)*T275*(-(y(1)*y(6)))/(y(1)*y(19)*y(1)*y(19))));
  g1(12,21)=1;
  g1(12,22)=(-(T138*(-params(8))));
  g1(12,23)=(-T460);
  g1(12,6)=(-((-y(22))*params(8)*T275*1/(y(1)*y(19))));
  g1(12,24)=(-(T405*(-(T124*T126))));
  g1(13,10)=(-(params(1)*(1-params(8))*(T149*y(43)*getPowerDeriv(y(10),(-params(4)),1)+T145*y(41)/y(24)*T314)));
  g1(13,19)=(-(T126*T124*T405));
  g1(13,41)=(-(params(1)*(1-params(8))*T145*T314*y(10)/y(24)));
  g1(13,22)=1;
  g1(13,43)=(-(params(1)*(1-params(8))*y(10)^(-params(4))*T149));
  g1(13,23)=(-(T124*T129*getPowerDeriv(y(23),1+params(3),1)+T126*T124*T460));
  g1(13,24)=(-(params(1)*(1-params(8))*T145*T314*(-(y(10)*y(41)))/(y(24)*y(24))+T126*T124*T405*(-(T124*T126))));
  g1(14,18)=(-(T165*1/y(23)));
  g1(14,19)=T160*T405;
  g1(14,21)=T129*y(24)*params(6)*y(23)^params(3)*(-1)/(y(21)*y(21));
  g1(14,23)=T160*T460+T129*1/y(21)*y(24)*params(6)*getPowerDeriv(y(23),params(3),1)-T165*(-y(18))/(y(23)*y(23));
  g1(14,24)=T160*T405*(-(T124*T126))+T129*1/y(21)*params(6)*y(23)^params(3);
  g1(15,1)=(-(y(19)^params(8)*(-y(6))/(y(1)*y(1))*T281));
  g1(15,19)=(-(T171*getPowerDeriv(y(19),params(8),1)));
  g1(15,6)=(-(y(19)^params(8)*T281*1/y(1)));
  g1(15,24)=1;
  g1(16,1)=(-((1-params(5))*(-y(7))/(y(1)*y(1))));
  g1(16,7)=(-((1-params(5))*1/y(1)));
  g1(16,25)=1;
  g1(16,27)=(-1);
  g1(17,1)=(-(log(y(35))*y(27)/y(8)));
  g1(17,10)=log(y(48))*y(42)*y(45)/y(27)*2*y(10)*y(45)/y(27);
  g1(17,42)=T191*log(y(48));
  g1(17,26)=1;
  g1(17,8)=(-(log(y(35))*T523));
  g1(17,27)=(-(log(y(35))*y(1)/y(8)-log(y(48))*y(42)*2*y(10)*y(45)/y(27)*(-(y(10)*y(45)))/(y(27)*y(27))));
  g1(17,45)=log(y(48))*y(42)*2*y(10)*y(45)/y(27)*y(10)/y(27);
  g1(17,34)=(-(1/y(34)));
  g1(17,35)=(-(y(1)*y(27)/y(8)*1/y(35)));
  g1(17,48)=y(42)*T191*1/y(48);
  g1(18,10)=(-(y(42)*params(2)*(params(7)-1)*y(40)/(params(7)*params(15)*y(25))));
  g1(18,40)=(-(y(42)*params(2)*y(10)*(params(7)-1)/(params(7)*params(15)*y(25))));
  g1(18,42)=(-(params(2)*y(10)*(params(7)-1)*y(40)/(params(7)*params(15)*y(25))+(1-params(5))*y(44)));
  g1(18,25)=(-(y(42)*(-(params(2)*y(10)*(params(7)-1)*y(40)*params(7)*params(15)))/(params(7)*params(15)*y(25)*params(7)*params(15)*y(25))));
  g1(18,26)=1;
  g1(18,44)=(-(y(42)*(1-params(5))));
  g1(19,9)=(-(params(12)*1/y(9)));
  g1(19,28)=1/y(28);
  g1(19,49)=(-params(13));
  g1(20,11)=(-y(13));
  g1(20,12)=(-y(13));
  g1(20,13)=(-(y(11)+y(12)-1));
  g1(20,14)=(-1);
  g1(20,25)=(-y(26));
  g1(20,26)=(-y(25));
  g1(20,29)=1;
  g1(20,30)=(-1);
  g1(21,10)=(-(y(42)*(y(37)*y(36)+y(46))));
  g1(21,36)=(-(y(37)*y(10)*y(42)));
  g1(21,37)=(-(y(10)*y(42)*y(36)));
  g1(21,42)=(-(y(10)*(y(37)*y(36)+y(46))));
  g1(21,30)=1;
  g1(21,46)=(-(y(10)*y(42)));
  g1(22,16)=(-1);
  g1(22,31)=1;
  g1(23,1)=(-(T243*T238*y(16)/y(4)*T290));
  g1(23,4)=(-(T243*T238*T290*T360));
  g1(23,16)=(-(T243*T238*T290*y(1)/y(4)));
  g1(23,32)=1;
  g1(24,1)=(-(T246*params(16)*y(16)/y(4)));
  g1(24,4)=(-(T246*params(16)*T360));
  g1(24,16)=(-(T246*params(16)*y(1)/y(4)));
  g1(24,33)=1;
  g1(25,1)=(-(T253*T249*y(27)/y(8)*T298));
  g1(25,8)=(-(T253*T249*T298*T523));
  g1(25,27)=(-(T253*T249*T298*y(1)/y(8)));
  g1(25,34)=1;
  g1(26,1)=(-(T256*params(17)*y(27)/y(8)));
  g1(26,8)=(-(T256*params(17)*T523));
  g1(26,27)=(-(T256*params(17)*y(1)/y(8)));
  g1(26,35)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],26,2401);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],26,117649);
end
end
