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
T57 = 1/params(8)*1/params(14);
T64 = (y(11)/y(16))^(1-params(11));
T79 = (y(10)*y(39)/y(16))^2;
T88 = (y(7)/y(1))^params(2);
T91 = y(23)^(1-params(2));
T110 = y(1)^(-params(4));
T117 = params(7)/(1+params(3));
T119 = y(23)^(1+params(3));
T121 = y(19)-y(24)*T117*T119;
T122 = T121^(-params(4));
T130 = (y(6)/(y(1)*y(19)))^(1-params(9));
T137 = y(10)^(-params(4))*y(43);
T141 = (y(10)*y(41)/y(24))^params(9);
T152 = y(24)*params(7)*y(23)^params(3)*1/y(21);
T157 = (1-params(2))*1/params(14)*(params(8)-1)/params(8);
T163 = (y(6)/y(1))^(1-params(9));
T182 = (y(10)*y(45)/y(27))^2;
T229 = params(15)/2;
T238 = params(16)/2;
T251 = getPowerDeriv(y(7)/y(1),params(2),1);
T262 = getPowerDeriv(y(6)/(y(1)*y(19)),1-params(9),1);
T268 = getPowerDeriv(y(6)/y(1),1-params(9),1);
T277 = 2*(y(1)*y(16)/y(4)-params(17));
T283 = 2*(y(1)*y(27)/y(8)-params(17));
T297 = getPowerDeriv(y(10)*y(41)/y(24),params(9),1);
T325 = getPowerDeriv(y(11)/y(16),1-params(11),1);
T344 = (-(y(1)*y(16)))/(y(4)*y(4));
T385 = getPowerDeriv(T121,(-params(4)),1);
T439 = T385*(-(y(24)*T117*getPowerDeriv(y(23),1+params(3),1)));
T502 = (-(y(1)*y(27)))/(y(8)*y(8));
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
rhs =1+y(33)*y(1)*y(16)/y(4)+y(32)-y(42)*y(47)*T79;
residual(7)= lhs-rhs;
lhs =y(18);
rhs =y(17);
residual(8)= lhs-rhs;
lhs =y(17);
rhs =T88*T91;
residual(9)= lhs-rhs;
lhs =y(18);
rhs =y(16)+y(19)+(1+y(34))*y(27);
residual(10)= lhs-rhs;
lhs =y(20);
rhs =params(1)*y(21)/y(5)*T110;
residual(11)= lhs-rhs;
lhs =y(21);
rhs =T122+y(22)*params(9)*T130;
residual(12)= lhs-rhs;
lhs =y(22);
rhs =params(1)*(1-params(9))*T137*T141-T119*T117*T122;
residual(13)= lhs-rhs;
lhs =T122*T152;
rhs =T157*y(18)/y(23);
residual(14)= lhs-rhs;
lhs =y(24);
rhs =y(19)^params(9)*T163;
residual(15)= lhs-rhs;
lhs =y(25);
rhs =y(27)+y(7)/y(1)*(1-params(6));
residual(16)= lhs-rhs;
lhs =y(26);
rhs =1+y(34)+y(1)*y(27)/y(8)*y(35)-y(42)*T182*y(48);
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
rhs =T229*(y(1)*y(16)/y(4)-params(17))^2;
residual(23)= lhs-rhs;
lhs =y(33);
rhs =params(15)*(y(1)*y(16)/y(4)-params(17));
residual(24)= lhs-rhs;
lhs =y(34);
rhs =T238*(y(1)*y(27)/y(8)-params(17))^2;
residual(25)= lhs-rhs;
lhs =y(35);
rhs =params(16)*(y(1)*y(27)/y(8)-params(17));
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
  g1(7,1)=(-(y(33)*y(16)/y(4)));
  g1(7,10)=y(42)*y(47)*y(39)/y(16)*2*y(10)*y(39)/y(16);
  g1(7,11)=y(28)*params(11)*y(13)*1/y(16)*T325;
  g1(7,13)=T64*y(28)*params(11);
  g1(7,4)=(-(y(33)*T344));
  g1(7,16)=y(28)*params(11)*y(13)*T325*(-y(11))/(y(16)*y(16))-(y(33)*y(1)/y(4)-y(42)*y(47)*2*y(10)*y(39)/y(16)*(-(y(10)*y(39)))/(y(16)*y(16)));
  g1(7,39)=y(42)*y(47)*2*y(10)*y(39)/y(16)*y(10)/y(16);
  g1(7,42)=y(47)*T79;
  g1(7,28)=params(11)*y(13)*T64;
  g1(7,32)=(-1);
  g1(7,33)=(-(y(1)*y(16)/y(4)));
  g1(7,47)=y(42)*T79;
  g1(8,17)=(-1);
  g1(8,18)=1;
  g1(9,1)=(-(T91*(-y(7))/(y(1)*y(1))*T251));
  g1(9,17)=1;
  g1(9,23)=(-(T88*getPowerDeriv(y(23),1-params(2),1)));
  g1(9,7)=(-(T91*T251*1/y(1)));
  g1(10,16)=(-1);
  g1(10,18)=1;
  g1(10,19)=(-1);
  g1(10,27)=(-(1+y(34)));
  g1(10,34)=(-y(27));
  g1(11,1)=(-(params(1)*y(21)/y(5)*getPowerDeriv(y(1),(-params(4)),1)));
  g1(11,20)=1;
  g1(11,5)=(-(T110*(-(params(1)*y(21)))/(y(5)*y(5))));
  g1(11,21)=(-(T110*params(1)/y(5)));
  g1(12,1)=(-(y(22)*params(9)*(-(y(19)*y(6)))/(y(1)*y(19)*y(1)*y(19))*T262));
  g1(12,19)=(-(T385+y(22)*params(9)*T262*(-(y(1)*y(6)))/(y(1)*y(19)*y(1)*y(19))));
  g1(12,21)=1;
  g1(12,22)=(-(params(9)*T130));
  g1(12,23)=(-T439);
  g1(12,6)=(-(y(22)*params(9)*T262*1/(y(1)*y(19))));
  g1(12,24)=(-(T385*(-(T117*T119))));
  g1(13,10)=(-(params(1)*(1-params(9))*(T141*y(43)*getPowerDeriv(y(10),(-params(4)),1)+T137*y(41)/y(24)*T297)));
  g1(13,19)=T119*T117*T385;
  g1(13,41)=(-(params(1)*(1-params(9))*T137*T297*y(10)/y(24)));
  g1(13,22)=1;
  g1(13,43)=(-(params(1)*(1-params(9))*y(10)^(-params(4))*T141));
  g1(13,23)=T117*T122*getPowerDeriv(y(23),1+params(3),1)+T119*T117*T439;
  g1(13,24)=(-(params(1)*(1-params(9))*T137*T297*(-(y(10)*y(41)))/(y(24)*y(24))-T119*T117*T385*(-(T117*T119))));
  g1(14,18)=(-(T157*1/y(23)));
  g1(14,19)=T152*T385;
  g1(14,21)=T122*y(24)*params(7)*y(23)^params(3)*(-1)/(y(21)*y(21));
  g1(14,23)=T152*T439+T122*1/y(21)*y(24)*params(7)*getPowerDeriv(y(23),params(3),1)-T157*(-y(18))/(y(23)*y(23));
  g1(14,24)=T152*T385*(-(T117*T119))+T122*1/y(21)*params(7)*y(23)^params(3);
  g1(15,1)=(-(y(19)^params(9)*(-y(6))/(y(1)*y(1))*T268));
  g1(15,19)=(-(T163*getPowerDeriv(y(19),params(9),1)));
  g1(15,6)=(-(y(19)^params(9)*T268*1/y(1)));
  g1(15,24)=1;
  g1(16,1)=(-((1-params(6))*(-y(7))/(y(1)*y(1))));
  g1(16,7)=(-((1-params(6))*1/y(1)));
  g1(16,25)=1;
  g1(16,27)=(-1);
  g1(17,1)=(-(y(35)*y(27)/y(8)));
  g1(17,10)=y(48)*y(42)*y(45)/y(27)*2*y(10)*y(45)/y(27);
  g1(17,42)=T182*y(48);
  g1(17,26)=1;
  g1(17,8)=(-(y(35)*T502));
  g1(17,27)=(-(y(35)*y(1)/y(8)-y(48)*y(42)*2*y(10)*y(45)/y(27)*(-(y(10)*y(45)))/(y(27)*y(27))));
  g1(17,45)=y(48)*y(42)*2*y(10)*y(45)/y(27)*y(10)/y(27);
  g1(17,34)=(-1);
  g1(17,35)=(-(y(1)*y(27)/y(8)));
  g1(17,48)=y(42)*T182;
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
  g1(23,1)=(-(T229*y(16)/y(4)*T277));
  g1(23,4)=(-(T229*T277*T344));
  g1(23,16)=(-(T229*T277*y(1)/y(4)));
  g1(23,32)=1;
  g1(24,1)=(-(params(15)*y(16)/y(4)));
  g1(24,4)=(-(params(15)*T344));
  g1(24,16)=(-(params(15)*y(1)/y(4)));
  g1(24,33)=1;
  g1(25,1)=(-(T238*y(27)/y(8)*T283));
  g1(25,8)=(-(T238*T283*T502));
  g1(25,27)=(-(T238*T283*y(1)/y(8)));
  g1(25,34)=1;
  g1(26,1)=(-(params(16)*y(27)/y(8)));
  g1(26,8)=(-(params(16)*T502));
  g1(26,27)=(-(params(16)*y(1)/y(8)));
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
