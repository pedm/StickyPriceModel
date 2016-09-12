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
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(27, 1);
T59 = 1/params(7)*1/params(16);
T67 = (y(12)/y(17))^(1-params(10));
T85 = (y(11)*y(41)/y(17))^2;
T94 = (y(7)/y(1))^params(2);
T97 = y(24)^(1-params(2));
T117 = y(1)^(-params(4));
T124 = params(6)/(1+params(3));
T126 = y(24)^(1+params(3));
T128 = y(20)-y(25)*T124*T126;
T129 = T128^(-params(4));
T138 = (y(6)/(y(1)*y(20)))^(1-params(8));
T145 = y(11)^(-params(4))*y(45);
T149 = (y(11)*y(43)/y(25))^params(8);
T160 = y(25)*params(6)*y(24)^params(3)*1/y(22);
T165 = (1-params(2))*1/params(16)*(params(7)-1)/params(7);
T171 = (y(6)/y(1))^(1-params(8));
T191 = (y(11)*y(47)/y(28))^2;
T244 = params(17)/2;
T249 = exp(T244*(y(1)*y(17)/y(4)-params(19))^2);
T252 = exp(params(17)*(y(1)*y(17)/y(4)-params(19)));
T255 = params(18)/2;
T259 = exp(T255*(y(1)*y(28)/y(8)-params(19))^2);
T262 = exp(params(18)*(y(1)*y(28)/y(8)-params(19)));
T278 = getPowerDeriv(y(7)/y(1),params(2),1);
T289 = getPowerDeriv(y(6)/(y(1)*y(20)),1-params(8),1);
T295 = getPowerDeriv(y(6)/y(1),1-params(8),1);
T304 = 2*(y(1)*y(17)/y(4)-params(19));
T312 = 2*(y(1)*y(28)/y(8)-params(19));
T328 = getPowerDeriv(y(11)*y(43)/y(25),params(8),1);
T356 = getPowerDeriv(y(12)/y(17),1-params(10),1);
T374 = (-(y(1)*y(17)))/(y(4)*y(4));
T419 = getPowerDeriv(T128,(-params(4)),1);
T474 = T419*(-(y(25)*T124*getPowerDeriv(y(24),1+params(3),1)));
T537 = (-(y(1)*y(28)))/(y(8)*y(8));
lhs =y(11);
rhs =params(9)*(1+params(11)*(y(12)+y(13)-1));
residual(1)= lhs-rhs;
lhs =y(12)*y(1);
rhs =params(9)*(y(2)+(1-params(11))*y(3));
residual(2)= lhs-rhs;
lhs =y(13);
rhs =params(15)*y(29)*y(12)^(1-params(10))*y(17)^params(10);
residual(3)= lhs-rhs;
lhs =y(14);
rhs =params(11)*y(15)+params(9)*(1-params(11))*y(44)*y(39);
residual(4)= lhs-rhs;
lhs =y(15);
rhs =y(16)+params(9)*y(44)*y(40);
residual(5)= lhs-rhs;
lhs =y(16);
rhs =T59*y(18);
residual(6)= lhs-rhs;
lhs =y(29)*y(14)*params(15)*params(10)*T67;
rhs =1+log(y(34))*y(1)*y(17)/y(4)+log(y(33))-y(44)*log(y(49))*T85;
residual(7)= lhs-rhs;
lhs =y(19);
rhs =y(18);
residual(8)= lhs-rhs;
lhs =y(18);
rhs =T94*T97;
residual(9)= lhs-rhs;
lhs =y(19);
rhs =y(17)+y(20)+(1+log(y(35)))*y(28);
residual(10)= lhs-rhs;
lhs =y(21);
rhs =params(1)*y(22)/y(5)*T117;
residual(11)= lhs-rhs;
lhs =y(22);
rhs =T129+(-y(23))*params(8)*T138;
residual(12)= lhs-rhs;
lhs =y(23);
rhs =params(1)*(1-params(8))*T145*T149+T126*T124*T129;
residual(13)= lhs-rhs;
lhs =T129*T160;
rhs =T165*y(19)/y(24);
residual(14)= lhs-rhs;
lhs =y(25);
rhs =y(20)^params(8)*T171;
residual(15)= lhs-rhs;
lhs =y(26);
rhs =y(28)+y(7)/y(1)*(1-params(5));
residual(16)= lhs-rhs;
lhs =y(27);
rhs =1+log(y(35))+y(1)*y(28)/y(8)*log(y(36))-y(44)*T191*log(y(50));
residual(17)= lhs-rhs;
lhs =y(27);
rhs =y(44)*(params(2)*y(11)*(params(7)-1)*y(42)/(params(7)*params(16)*y(26))+(1-params(5))*y(46));
residual(18)= lhs-rhs;
lhs =log(y(29));
rhs =params(14)*x(it_, 1)+params(12)*log(y(9))-params(13)*log(y(10));
residual(19)= lhs-rhs;
lhs =y(30);
rhs =y(15)+y(26)*y(27)+(y(12)+y(13)-1)*y(14)+y(31);
residual(20)= lhs-rhs;
lhs =y(31);
rhs =y(11)*y(44)*(y(39)*y(38)+y(48));
residual(21)= lhs-rhs;
lhs =y(32);
rhs =y(17);
residual(22)= lhs-rhs;
lhs =y(33);
rhs =T249;
residual(23)= lhs-rhs;
lhs =y(34);
rhs =T252;
residual(24)= lhs-rhs;
lhs =y(35);
rhs =T259;
residual(25)= lhs-rhs;
lhs =y(36);
rhs =T262;
residual(26)= lhs-rhs;
lhs =y(37);
rhs =y(9);
residual(27)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(27, 51);

  %
  % Jacobian matrix
  %

  g1(1,11)=1;
  g1(1,12)=(-(params(9)*params(11)));
  g1(1,13)=(-(params(9)*params(11)));
  g1(2,1)=y(12);
  g1(2,2)=(-params(9));
  g1(2,12)=y(1);
  g1(2,3)=(-(params(9)*(1-params(11))));
  g1(3,12)=(-(y(17)^params(10)*params(15)*y(29)*getPowerDeriv(y(12),1-params(10),1)));
  g1(3,13)=1;
  g1(3,17)=(-(params(15)*y(29)*y(12)^(1-params(10))*getPowerDeriv(y(17),params(10),1)));
  g1(3,29)=(-(y(17)^params(10)*params(15)*y(12)^(1-params(10))));
  g1(4,14)=1;
  g1(4,39)=(-(params(9)*(1-params(11))*y(44)));
  g1(4,15)=(-params(11));
  g1(4,44)=(-(params(9)*(1-params(11))*y(39)));
  g1(5,15)=1;
  g1(5,40)=(-(params(9)*y(44)));
  g1(5,16)=(-1);
  g1(5,44)=(-(params(9)*y(40)));
  g1(6,16)=1;
  g1(6,18)=(-T59);
  g1(7,1)=(-(log(y(34))*y(17)/y(4)));
  g1(7,11)=y(44)*log(y(49))*y(41)/y(17)*2*y(11)*y(41)/y(17);
  g1(7,12)=y(29)*y(14)*params(15)*params(10)*1/y(17)*T356;
  g1(7,14)=T67*y(29)*params(15)*params(10);
  g1(7,4)=(-(log(y(34))*T374));
  g1(7,17)=y(29)*y(14)*params(15)*params(10)*T356*(-y(12))/(y(17)*y(17))-(log(y(34))*y(1)/y(4)-y(44)*log(y(49))*2*y(11)*y(41)/y(17)*(-(y(11)*y(41)))/(y(17)*y(17)));
  g1(7,41)=y(44)*log(y(49))*2*y(11)*y(41)/y(17)*y(11)/y(17);
  g1(7,44)=log(y(49))*T85;
  g1(7,29)=y(14)*params(15)*params(10)*T67;
  g1(7,33)=(-(1/y(33)));
  g1(7,34)=(-(y(1)*y(17)/y(4)*1/y(34)));
  g1(7,49)=T85*y(44)*1/y(49);
  g1(8,18)=(-1);
  g1(8,19)=1;
  g1(9,1)=(-(T97*(-y(7))/(y(1)*y(1))*T278));
  g1(9,18)=1;
  g1(9,24)=(-(T94*getPowerDeriv(y(24),1-params(2),1)));
  g1(9,7)=(-(T97*T278*1/y(1)));
  g1(10,17)=(-1);
  g1(10,19)=1;
  g1(10,20)=(-1);
  g1(10,28)=(-(1+log(y(35))));
  g1(10,35)=(-(y(28)*1/y(35)));
  g1(11,1)=(-(params(1)*y(22)/y(5)*getPowerDeriv(y(1),(-params(4)),1)));
  g1(11,21)=1;
  g1(11,5)=(-(T117*(-(params(1)*y(22)))/(y(5)*y(5))));
  g1(11,22)=(-(T117*params(1)/y(5)));
  g1(12,1)=(-((-y(23))*params(8)*(-(y(20)*y(6)))/(y(1)*y(20)*y(1)*y(20))*T289));
  g1(12,20)=(-(T419+(-y(23))*params(8)*T289*(-(y(1)*y(6)))/(y(1)*y(20)*y(1)*y(20))));
  g1(12,22)=1;
  g1(12,23)=(-(T138*(-params(8))));
  g1(12,24)=(-T474);
  g1(12,6)=(-((-y(23))*params(8)*T289*1/(y(1)*y(20))));
  g1(12,25)=(-(T419*(-(T124*T126))));
  g1(13,11)=(-(params(1)*(1-params(8))*(T149*y(45)*getPowerDeriv(y(11),(-params(4)),1)+T145*y(43)/y(25)*T328)));
  g1(13,20)=(-(T126*T124*T419));
  g1(13,43)=(-(params(1)*(1-params(8))*T145*T328*y(11)/y(25)));
  g1(13,23)=1;
  g1(13,45)=(-(params(1)*(1-params(8))*y(11)^(-params(4))*T149));
  g1(13,24)=(-(T124*T129*getPowerDeriv(y(24),1+params(3),1)+T126*T124*T474));
  g1(13,25)=(-(params(1)*(1-params(8))*T145*T328*(-(y(11)*y(43)))/(y(25)*y(25))+T126*T124*T419*(-(T124*T126))));
  g1(14,19)=(-(T165*1/y(24)));
  g1(14,20)=T160*T419;
  g1(14,22)=T129*y(25)*params(6)*y(24)^params(3)*(-1)/(y(22)*y(22));
  g1(14,24)=T160*T474+T129*1/y(22)*y(25)*params(6)*getPowerDeriv(y(24),params(3),1)-T165*(-y(19))/(y(24)*y(24));
  g1(14,25)=T160*T419*(-(T124*T126))+T129*1/y(22)*params(6)*y(24)^params(3);
  g1(15,1)=(-(y(20)^params(8)*(-y(6))/(y(1)*y(1))*T295));
  g1(15,20)=(-(T171*getPowerDeriv(y(20),params(8),1)));
  g1(15,6)=(-(y(20)^params(8)*T295*1/y(1)));
  g1(15,25)=1;
  g1(16,1)=(-((1-params(5))*(-y(7))/(y(1)*y(1))));
  g1(16,7)=(-((1-params(5))*1/y(1)));
  g1(16,26)=1;
  g1(16,28)=(-1);
  g1(17,1)=(-(log(y(36))*y(28)/y(8)));
  g1(17,11)=log(y(50))*y(44)*y(47)/y(28)*2*y(11)*y(47)/y(28);
  g1(17,44)=T191*log(y(50));
  g1(17,27)=1;
  g1(17,8)=(-(log(y(36))*T537));
  g1(17,28)=(-(log(y(36))*y(1)/y(8)-log(y(50))*y(44)*2*y(11)*y(47)/y(28)*(-(y(11)*y(47)))/(y(28)*y(28))));
  g1(17,47)=log(y(50))*y(44)*2*y(11)*y(47)/y(28)*y(11)/y(28);
  g1(17,35)=(-(1/y(35)));
  g1(17,36)=(-(y(1)*y(28)/y(8)*1/y(36)));
  g1(17,50)=y(44)*T191*1/y(50);
  g1(18,11)=(-(y(44)*params(2)*(params(7)-1)*y(42)/(params(7)*params(16)*y(26))));
  g1(18,42)=(-(y(44)*params(2)*y(11)*(params(7)-1)/(params(7)*params(16)*y(26))));
  g1(18,44)=(-(params(2)*y(11)*(params(7)-1)*y(42)/(params(7)*params(16)*y(26))+(1-params(5))*y(46)));
  g1(18,26)=(-(y(44)*(-(params(2)*y(11)*(params(7)-1)*y(42)*params(7)*params(16)))/(params(7)*params(16)*y(26)*params(7)*params(16)*y(26))));
  g1(18,27)=1;
  g1(18,46)=(-(y(44)*(1-params(5))));
  g1(19,9)=(-(params(12)*1/y(9)));
  g1(19,29)=1/y(29);
  g1(19,51)=(-params(14));
  g1(19,10)=params(13)*1/y(10);
  g1(20,12)=(-y(14));
  g1(20,13)=(-y(14));
  g1(20,14)=(-(y(12)+y(13)-1));
  g1(20,15)=(-1);
  g1(20,26)=(-y(27));
  g1(20,27)=(-y(26));
  g1(20,30)=1;
  g1(20,31)=(-1);
  g1(21,11)=(-(y(44)*(y(39)*y(38)+y(48))));
  g1(21,38)=(-(y(39)*y(11)*y(44)));
  g1(21,39)=(-(y(11)*y(44)*y(38)));
  g1(21,44)=(-(y(11)*(y(39)*y(38)+y(48))));
  g1(21,31)=1;
  g1(21,48)=(-(y(11)*y(44)));
  g1(22,17)=(-1);
  g1(22,32)=1;
  g1(23,1)=(-(T249*T244*y(17)/y(4)*T304));
  g1(23,4)=(-(T249*T244*T304*T374));
  g1(23,17)=(-(T249*T244*T304*y(1)/y(4)));
  g1(23,33)=1;
  g1(24,1)=(-(T252*params(17)*y(17)/y(4)));
  g1(24,4)=(-(T252*params(17)*T374));
  g1(24,17)=(-(T252*params(17)*y(1)/y(4)));
  g1(24,34)=1;
  g1(25,1)=(-(T259*T255*y(28)/y(8)*T312));
  g1(25,8)=(-(T259*T255*T312*T537));
  g1(25,28)=(-(T259*T255*T312*y(1)/y(8)));
  g1(25,35)=1;
  g1(26,1)=(-(T262*params(18)*y(28)/y(8)));
  g1(26,8)=(-(T262*params(18)*T537));
  g1(26,28)=(-(T262*params(18)*y(1)/y(8)));
  g1(26,36)=1;
  g1(27,9)=(-1);
  g1(27,37)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],27,2601);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],27,132651);
end
end
