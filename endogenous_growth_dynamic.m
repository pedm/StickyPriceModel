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

residual = zeros(26, 1);
T57 = 1/params(8)*1/params(14);
T64 = (y(11)/y(16))^(1-params(11));
T82 = (y(10)*y(39)/y(16))^2;
T91 = (y(7)/y(1))^params(2);
T94 = y(23)^(1-params(2));
T114 = y(1)^(-params(4));
T121 = params(7)/(1+params(3));
T123 = y(23)^(1+params(3));
T125 = y(19)-y(24)*T121*T123;
T126 = T125^(-params(4));
T135 = (y(6)/(y(1)*y(19)))^(1-params(9));
T142 = y(10)^(-params(4))*y(43);
T146 = (y(10)*y(41)/y(24))^params(9);
T157 = y(24)*params(7)*y(23)^params(3)*1/y(21);
T162 = (1-params(2))*1/params(14)*(params(8)-1)/params(8);
T168 = (y(6)/y(1))^(1-params(9));
T188 = (y(10)*y(45)/y(27))^2;
T236 = params(15)/2;
T241 = exp(T236*(y(1)*y(16)/y(4)-params(17))^2);
T244 = exp(params(15)*(y(1)*y(16)/y(4)-params(17)));
T247 = params(16)/2;
T251 = exp(T247*(y(1)*y(27)/y(8)-params(17))^2);
T254 = exp(params(16)*(y(1)*y(27)/y(8)-params(17)));
T262 = getPowerDeriv(y(7)/y(1),params(2),1);
T273 = getPowerDeriv(y(6)/(y(1)*y(19)),1-params(9),1);
T279 = getPowerDeriv(y(6)/y(1),1-params(9),1);
T288 = 2*(y(1)*y(16)/y(4)-params(17));
T296 = 2*(y(1)*y(27)/y(8)-params(17));
T312 = getPowerDeriv(y(10)*y(41)/y(24),params(9),1);
T340 = getPowerDeriv(y(11)/y(16),1-params(11),1);
T359 = (-(y(1)*y(16)))/(y(4)*y(4));
T404 = getPowerDeriv(T125,(-params(4)),1);
T459 = T404*(-(y(24)*T121*getPowerDeriv(y(23),1+params(3),1)));
T522 = (-(y(1)*y(27)))/(y(8)*y(8));
lhs =y(10);
rhs =params(10)*(1+params(12)*(y(11)+y(12)-1));
residual(1)= lhs-rhs;
lhs =y(11)*y(1);
rhs =params(10)*(y(2)+(1-params(12))*y(3));
residual(2)= lhs-rhs;
lhs =y(12);
rhs =y(28)*y(11)^(1-params(11))*y(16)^params(11);
residual(3)= lhs-rhs;
lhs =y(13);
rhs =params(12)*y(14)+params(10)*(1-params(12))*y(42)*y(37);
residual(4)= lhs-rhs;
lhs =y(14);
rhs =y(15)+params(10)*y(42)*y(38);
residual(5)= lhs-rhs;
lhs =y(15);
rhs =T57*y(17);
residual(6)= lhs-rhs;
lhs =y(28)*params(11)*y(13)*T64;
rhs =1+log(y(33))*y(1)*y(16)/y(4)+log(y(32))-y(42)*log(y(47))*T82;
residual(7)= lhs-rhs;
lhs =y(18);
rhs =y(17);
residual(8)= lhs-rhs;
lhs =y(17);
rhs =T91*T94;
residual(9)= lhs-rhs;
lhs =y(18);
rhs =y(16)+y(19)+(1+log(y(34)))*y(27);
residual(10)= lhs-rhs;
lhs =y(20);
rhs =params(1)*y(21)/y(5)*T114;
residual(11)= lhs-rhs;
lhs =y(21);
rhs =T126+(-y(22))*params(9)*T135;
residual(12)= lhs-rhs;
lhs =y(22);
rhs =params(1)*(1-params(9))*T142*T146+T123*T121*T126;
residual(13)= lhs-rhs;
lhs =T126*T157;
rhs =T162*y(18)/y(23);
residual(14)= lhs-rhs;
lhs =y(24);
rhs =y(19)^params(9)*T168;
residual(15)= lhs-rhs;
lhs =y(25);
rhs =y(27)+y(7)/y(1)*(1-params(6));
residual(16)= lhs-rhs;
lhs =y(26);
rhs =1+log(y(34))+y(1)*y(27)/y(8)*log(y(35))-y(42)*T188*log(y(48));
residual(17)= lhs-rhs;
lhs =y(26);
rhs =y(42)*(params(2)*y(10)*(params(8)-1)*y(40)/(params(8)*params(14)*y(25))+(1-params(6))*y(44));
residual(18)= lhs-rhs;
lhs =log(y(28));
rhs =params(13)*log(y(9))+0.1*x(it_, 1);
residual(19)= lhs-rhs;
lhs =y(29);
rhs =y(14)+y(25)*y(26)+(y(11)+y(12)-1)*y(13)+y(30);
residual(20)= lhs-rhs;
lhs =y(30);
rhs =y(10)*y(42)*(y(37)*y(36)+y(46));
residual(21)= lhs-rhs;
lhs =y(31);
rhs =y(12)*y(13);
residual(22)= lhs-rhs;
lhs =y(32);
rhs =T241;
residual(23)= lhs-rhs;
lhs =y(33);
rhs =T244;
residual(24)= lhs-rhs;
lhs =y(34);
rhs =T251;
residual(25)= lhs-rhs;
lhs =y(35);
rhs =T254;
residual(26)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(26, 49);

  %
  % Jacobian matrix
  %

  g1(1,10)=1;
  g1(1,11)=(-(params(10)*params(12)));
  g1(1,12)=(-(params(10)*params(12)));
  g1(2,1)=y(11);
  g1(2,2)=(-params(10));
  g1(2,11)=y(1);
  g1(2,3)=(-(params(10)*(1-params(12))));
  g1(3,11)=(-(y(16)^params(11)*y(28)*getPowerDeriv(y(11),1-params(11),1)));
  g1(3,12)=1;
  g1(3,16)=(-(y(28)*y(11)^(1-params(11))*getPowerDeriv(y(16),params(11),1)));
  g1(3,28)=(-(y(11)^(1-params(11))*y(16)^params(11)));
  g1(4,13)=1;
  g1(4,37)=(-(params(10)*(1-params(12))*y(42)));
  g1(4,14)=(-params(12));
  g1(4,42)=(-(params(10)*(1-params(12))*y(37)));
  g1(5,14)=1;
  g1(5,38)=(-(params(10)*y(42)));
  g1(5,15)=(-1);
  g1(5,42)=(-(params(10)*y(38)));
  g1(6,15)=1;
  g1(6,17)=(-T57);
  g1(7,1)=(-(log(y(33))*y(16)/y(4)));
  g1(7,10)=y(42)*log(y(47))*y(39)/y(16)*2*y(10)*y(39)/y(16);
  g1(7,11)=y(28)*params(11)*y(13)*1/y(16)*T340;
  g1(7,13)=T64*y(28)*params(11);
  g1(7,4)=(-(log(y(33))*T359));
  g1(7,16)=y(28)*params(11)*y(13)*T340*(-y(11))/(y(16)*y(16))-(log(y(33))*y(1)/y(4)-y(42)*log(y(47))*2*y(10)*y(39)/y(16)*(-(y(10)*y(39)))/(y(16)*y(16)));
  g1(7,39)=y(42)*log(y(47))*2*y(10)*y(39)/y(16)*y(10)/y(16);
  g1(7,42)=log(y(47))*T82;
  g1(7,28)=params(11)*y(13)*T64;
  g1(7,32)=(-(1/y(32)));
  g1(7,33)=(-(y(1)*y(16)/y(4)*1/y(33)));
  g1(7,47)=T82*y(42)*1/y(47);
  g1(8,17)=(-1);
  g1(8,18)=1;
  g1(9,1)=(-(T94*(-y(7))/(y(1)*y(1))*T262));
  g1(9,17)=1;
  g1(9,23)=(-(T91*getPowerDeriv(y(23),1-params(2),1)));
  g1(9,7)=(-(T94*T262*1/y(1)));
  g1(10,16)=(-1);
  g1(10,18)=1;
  g1(10,19)=(-1);
  g1(10,27)=(-(1+log(y(34))));
  g1(10,34)=(-(y(27)*1/y(34)));
  g1(11,1)=(-(params(1)*y(21)/y(5)*getPowerDeriv(y(1),(-params(4)),1)));
  g1(11,20)=1;
  g1(11,5)=(-(T114*(-(params(1)*y(21)))/(y(5)*y(5))));
  g1(11,21)=(-(T114*params(1)/y(5)));
  g1(12,1)=(-((-y(22))*params(9)*(-(y(19)*y(6)))/(y(1)*y(19)*y(1)*y(19))*T273));
  g1(12,19)=(-(T404+(-y(22))*params(9)*T273*(-(y(1)*y(6)))/(y(1)*y(19)*y(1)*y(19))));
  g1(12,21)=1;
  g1(12,22)=(-(T135*(-params(9))));
  g1(12,23)=(-T459);
  g1(12,6)=(-((-y(22))*params(9)*T273*1/(y(1)*y(19))));
  g1(12,24)=(-(T404*(-(T121*T123))));
  g1(13,10)=(-(params(1)*(1-params(9))*(T146*y(43)*getPowerDeriv(y(10),(-params(4)),1)+T142*y(41)/y(24)*T312)));
  g1(13,19)=(-(T123*T121*T404));
  g1(13,41)=(-(params(1)*(1-params(9))*T142*T312*y(10)/y(24)));
  g1(13,22)=1;
  g1(13,43)=(-(params(1)*(1-params(9))*y(10)^(-params(4))*T146));
  g1(13,23)=(-(T121*T126*getPowerDeriv(y(23),1+params(3),1)+T123*T121*T459));
  g1(13,24)=(-(params(1)*(1-params(9))*T142*T312*(-(y(10)*y(41)))/(y(24)*y(24))+T123*T121*T404*(-(T121*T123))));
  g1(14,18)=(-(T162*1/y(23)));
  g1(14,19)=T157*T404;
  g1(14,21)=T126*y(24)*params(7)*y(23)^params(3)*(-1)/(y(21)*y(21));
  g1(14,23)=T157*T459+T126*1/y(21)*y(24)*params(7)*getPowerDeriv(y(23),params(3),1)-T162*(-y(18))/(y(23)*y(23));
  g1(14,24)=T157*T404*(-(T121*T123))+T126*1/y(21)*params(7)*y(23)^params(3);
  g1(15,1)=(-(y(19)^params(9)*(-y(6))/(y(1)*y(1))*T279));
  g1(15,19)=(-(T168*getPowerDeriv(y(19),params(9),1)));
  g1(15,6)=(-(y(19)^params(9)*T279*1/y(1)));
  g1(15,24)=1;
  g1(16,1)=(-((1-params(6))*(-y(7))/(y(1)*y(1))));
  g1(16,7)=(-((1-params(6))*1/y(1)));
  g1(16,25)=1;
  g1(16,27)=(-1);
  g1(17,1)=(-(log(y(35))*y(27)/y(8)));
  g1(17,10)=log(y(48))*y(42)*y(45)/y(27)*2*y(10)*y(45)/y(27);
  g1(17,42)=T188*log(y(48));
  g1(17,26)=1;
  g1(17,8)=(-(log(y(35))*T522));
  g1(17,27)=(-(log(y(35))*y(1)/y(8)-log(y(48))*y(42)*2*y(10)*y(45)/y(27)*(-(y(10)*y(45)))/(y(27)*y(27))));
  g1(17,45)=log(y(48))*y(42)*2*y(10)*y(45)/y(27)*y(10)/y(27);
  g1(17,34)=(-(1/y(34)));
  g1(17,35)=(-(y(1)*y(27)/y(8)*1/y(35)));
  g1(17,48)=y(42)*T188*1/y(48);
  g1(18,10)=(-(y(42)*params(2)*(params(8)-1)*y(40)/(params(8)*params(14)*y(25))));
  g1(18,40)=(-(y(42)*params(2)*y(10)*(params(8)-1)/(params(8)*params(14)*y(25))));
  g1(18,42)=(-(params(2)*y(10)*(params(8)-1)*y(40)/(params(8)*params(14)*y(25))+(1-params(6))*y(44)));
  g1(18,25)=(-(y(42)*(-(params(2)*y(10)*(params(8)-1)*y(40)*params(8)*params(14)))/(params(8)*params(14)*y(25)*params(8)*params(14)*y(25))));
  g1(18,26)=1;
  g1(18,44)=(-(y(42)*(1-params(6))));
  g1(19,9)=(-(params(13)*1/y(9)));
  g1(19,28)=1/y(28);
  g1(19,49)=(-0.1);
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
  g1(22,12)=(-y(13));
  g1(22,13)=(-y(12));
  g1(22,31)=1;
  g1(23,1)=(-(T241*T236*y(16)/y(4)*T288));
  g1(23,4)=(-(T241*T236*T288*T359));
  g1(23,16)=(-(T241*T236*T288*y(1)/y(4)));
  g1(23,32)=1;
  g1(24,1)=(-(T244*params(15)*y(16)/y(4)));
  g1(24,4)=(-(T244*params(15)*T359));
  g1(24,16)=(-(T244*params(15)*y(1)/y(4)));
  g1(24,33)=1;
  g1(25,1)=(-(T251*T247*y(27)/y(8)*T296));
  g1(25,8)=(-(T251*T247*T296*T522));
  g1(25,27)=(-(T251*T247*T296*y(1)/y(8)));
  g1(25,34)=1;
  g1(26,1)=(-(T254*params(16)*y(27)/y(8)));
  g1(26,8)=(-(T254*params(16)*T522));
  g1(26,27)=(-(T254*params(16)*y(1)/y(8)));
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
