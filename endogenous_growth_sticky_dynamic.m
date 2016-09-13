function [residual, g1, g2, g3] = endogenous_growth_sticky_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(33, 1);
T58 = 1/y(33);
T67 = (y(12)/y(17))^(1-params(10));
T85 = (y(11)*y(47)/y(17))^2;
T94 = (y(7)/y(1))^params(2);
T97 = y(24)^(1-params(2));
T117 = y(1)^(-params(4));
T124 = params(6)/(1+params(3));
T126 = y(24)^(1+params(3));
T128 = y(20)-y(25)*T124*T126;
T129 = T128^(-params(4));
T138 = (y(6)/(y(1)*y(20)))^(1-params(8));
T145 = y(11)^(-params(4))*y(51);
T149 = (y(11)*y(49)/y(25))^params(8);
T160 = y(25)*params(6)*y(24)^params(3)*1/y(22);
T165 = (1-params(2))*T58*(params(7)-1)/params(7);
T171 = (y(6)/y(1))^(1-params(8));
T191 = (y(11)*y(53)/y(28))^2;
T244 = params(17)/2;
T249 = exp(T244*(y(1)*y(17)/y(4)-params(19))^2);
T252 = exp(params(17)*(y(1)*y(17)/y(4)-params(19)));
T255 = params(18)/2;
T259 = exp(T255*(y(1)*y(28)/y(8)-params(19))^2);
T262 = exp(params(18)*(y(1)*y(28)/y(8)-params(19)));
T280 = params(23)/(params(23)-1)*y(36)/y(37);
T288 = params(1)*params(24)*y(11)^(1-params(4));
T313 = T58/(1/params(16));
T332 = getPowerDeriv(y(7)/y(1),params(2),1);
T343 = getPowerDeriv(y(6)/(y(1)*y(20)),1-params(8),1);
T349 = getPowerDeriv(y(6)/y(1),1-params(8),1);
T358 = 2*(y(1)*y(17)/y(4)-params(19));
T366 = 2*(y(1)*y(28)/y(8)-params(19));
T382 = getPowerDeriv(y(11)*y(49)/y(25),params(8),1);
T418 = getPowerDeriv(y(12)/y(17),1-params(10),1);
T436 = (-(y(1)*y(17)))/(y(4)*y(4));
T483 = getPowerDeriv(T128,(-params(4)),1);
T543 = T483*(-(y(25)*T124*getPowerDeriv(y(24),1+params(3),1)));
T605 = (-(y(1)*y(28)))/(y(8)*y(8));
T648 = (-1)/(y(33)*y(33));
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
rhs =params(11)*y(15)+params(9)*(1-params(11))*y(50)*y(45);
residual(4)= lhs-rhs;
lhs =y(15);
rhs =y(16)+params(9)*y(50)*y(46);
residual(5)= lhs-rhs;
lhs =y(16);
rhs =1/params(7)*T58*y(18);
residual(6)= lhs-rhs;
lhs =y(29)*y(14)*params(15)*params(10)*T67;
rhs =1+log(y(40))*y(1)*y(17)/y(4)+log(y(39))-y(50)*log(y(58))*T85;
residual(7)= lhs-rhs;
lhs =y(19);
rhs =y(18);
residual(8)= lhs-rhs;
lhs =y(18);
rhs =T94*T97;
residual(9)= lhs-rhs;
lhs =y(19);
rhs =y(17)+y(20)+(1+log(y(41)))*y(28);
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
rhs =1+log(y(41))+y(1)*y(28)/y(8)*log(y(42))-y(50)*T191*log(y(59));
residual(17)= lhs-rhs;
lhs =y(27);
rhs =y(50)*(params(2)*y(11)*(params(7)-1)*y(48)/(params(7)*y(33)*y(26))+(1-params(5))*y(52));
residual(18)= lhs-rhs;
lhs =log(y(29));
rhs =params(14)*x(it_, 1)+params(12)*log(y(9))-params(13)*log(y(10));
residual(19)= lhs-rhs;
lhs =y(30);
rhs =y(15)+y(26)*y(27)+(y(12)+y(13)-1)*y(14)+y(31);
residual(20)= lhs-rhs;
lhs =y(31);
rhs =y(11)*y(50)*(y(45)*y(44)+y(54));
residual(21)= lhs-rhs;
lhs =y(32);
rhs =y(17);
residual(22)= lhs-rhs;
lhs =y(39);
rhs =T249;
residual(23)= lhs-rhs;
lhs =y(40);
rhs =T252;
residual(24)= lhs-rhs;
lhs =y(41);
rhs =T259;
residual(25)= lhs-rhs;
lhs =y(42);
rhs =T262;
residual(26)= lhs-rhs;
lhs =y(34)^(1-params(23));
rhs =params(24)+(1-params(24))*y(35)^(1-params(23));
residual(27)= lhs-rhs;
lhs =y(35);
rhs =y(34)*T280;
residual(28)= lhs-rhs;
lhs =y(36);
rhs =y(19)*T58*y(22)+T288*y(55)^params(23)*y(56);
residual(29)= lhs-rhs;
lhs =y(37);
rhs =y(19)*y(22)+T288*y(55)^(params(23)-1)*y(57);
residual(30)= lhs-rhs;
lhs =1;
rhs =y(50)*y(38)/y(55);
residual(31)= lhs-rhs;
lhs =y(38)/params(20);
rhs =y(34)^params(21)*T313^params(22);
residual(32)= lhs-rhs;
lhs =y(43);
rhs =y(9);
residual(33)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(33, 60);

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
  g1(4,45)=(-(params(9)*(1-params(11))*y(50)));
  g1(4,15)=(-params(11));
  g1(4,50)=(-(params(9)*(1-params(11))*y(45)));
  g1(5,15)=1;
  g1(5,46)=(-(params(9)*y(50)));
  g1(5,16)=(-1);
  g1(5,50)=(-(params(9)*y(46)));
  g1(6,16)=1;
  g1(6,18)=(-(1/params(7)*T58));
  g1(6,33)=(-(y(18)*1/params(7)*T648));
  g1(7,1)=(-(log(y(40))*y(17)/y(4)));
  g1(7,11)=y(50)*log(y(58))*y(47)/y(17)*2*y(11)*y(47)/y(17);
  g1(7,12)=y(29)*y(14)*params(15)*params(10)*1/y(17)*T418;
  g1(7,14)=T67*y(29)*params(15)*params(10);
  g1(7,4)=(-(log(y(40))*T436));
  g1(7,17)=y(29)*y(14)*params(15)*params(10)*T418*(-y(12))/(y(17)*y(17))-(log(y(40))*y(1)/y(4)-y(50)*log(y(58))*2*y(11)*y(47)/y(17)*(-(y(11)*y(47)))/(y(17)*y(17)));
  g1(7,47)=y(50)*log(y(58))*2*y(11)*y(47)/y(17)*y(11)/y(17);
  g1(7,50)=log(y(58))*T85;
  g1(7,29)=y(14)*params(15)*params(10)*T67;
  g1(7,39)=(-(1/y(39)));
  g1(7,40)=(-(y(1)*y(17)/y(4)*1/y(40)));
  g1(7,58)=T85*y(50)*1/y(58);
  g1(8,18)=(-1);
  g1(8,19)=1;
  g1(9,1)=(-(T97*(-y(7))/(y(1)*y(1))*T332));
  g1(9,18)=1;
  g1(9,24)=(-(T94*getPowerDeriv(y(24),1-params(2),1)));
  g1(9,7)=(-(T97*T332*1/y(1)));
  g1(10,17)=(-1);
  g1(10,19)=1;
  g1(10,20)=(-1);
  g1(10,28)=(-(1+log(y(41))));
  g1(10,41)=(-(y(28)*1/y(41)));
  g1(11,1)=(-(params(1)*y(22)/y(5)*getPowerDeriv(y(1),(-params(4)),1)));
  g1(11,21)=1;
  g1(11,5)=(-(T117*(-(params(1)*y(22)))/(y(5)*y(5))));
  g1(11,22)=(-(T117*params(1)/y(5)));
  g1(12,1)=(-((-y(23))*params(8)*(-(y(20)*y(6)))/(y(1)*y(20)*y(1)*y(20))*T343));
  g1(12,20)=(-(T483+(-y(23))*params(8)*T343*(-(y(1)*y(6)))/(y(1)*y(20)*y(1)*y(20))));
  g1(12,22)=1;
  g1(12,23)=(-(T138*(-params(8))));
  g1(12,24)=(-T543);
  g1(12,6)=(-((-y(23))*params(8)*T343*1/(y(1)*y(20))));
  g1(12,25)=(-(T483*(-(T124*T126))));
  g1(13,11)=(-(params(1)*(1-params(8))*(T149*y(51)*getPowerDeriv(y(11),(-params(4)),1)+T145*y(49)/y(25)*T382)));
  g1(13,20)=(-(T126*T124*T483));
  g1(13,49)=(-(params(1)*(1-params(8))*T145*T382*y(11)/y(25)));
  g1(13,23)=1;
  g1(13,51)=(-(params(1)*(1-params(8))*y(11)^(-params(4))*T149));
  g1(13,24)=(-(T124*T129*getPowerDeriv(y(24),1+params(3),1)+T126*T124*T543));
  g1(13,25)=(-(params(1)*(1-params(8))*T145*T382*(-(y(11)*y(49)))/(y(25)*y(25))+T126*T124*T483*(-(T124*T126))));
  g1(14,19)=(-(T165*1/y(24)));
  g1(14,20)=T160*T483;
  g1(14,22)=T129*y(25)*params(6)*y(24)^params(3)*(-1)/(y(22)*y(22));
  g1(14,24)=T160*T543+T129*1/y(22)*y(25)*params(6)*getPowerDeriv(y(24),params(3),1)-T165*(-y(19))/(y(24)*y(24));
  g1(14,25)=T160*T483*(-(T124*T126))+T129*1/y(22)*params(6)*y(24)^params(3);
  g1(14,33)=(-(y(19)/y(24)*(1-params(2))*(params(7)-1)/params(7)*T648));
  g1(15,1)=(-(y(20)^params(8)*(-y(6))/(y(1)*y(1))*T349));
  g1(15,20)=(-(T171*getPowerDeriv(y(20),params(8),1)));
  g1(15,6)=(-(y(20)^params(8)*T349*1/y(1)));
  g1(15,25)=1;
  g1(16,1)=(-((1-params(5))*(-y(7))/(y(1)*y(1))));
  g1(16,7)=(-((1-params(5))*1/y(1)));
  g1(16,26)=1;
  g1(16,28)=(-1);
  g1(17,1)=(-(log(y(42))*y(28)/y(8)));
  g1(17,11)=log(y(59))*y(50)*y(53)/y(28)*2*y(11)*y(53)/y(28);
  g1(17,50)=T191*log(y(59));
  g1(17,27)=1;
  g1(17,8)=(-(log(y(42))*T605));
  g1(17,28)=(-(log(y(42))*y(1)/y(8)-log(y(59))*y(50)*2*y(11)*y(53)/y(28)*(-(y(11)*y(53)))/(y(28)*y(28))));
  g1(17,53)=log(y(59))*y(50)*2*y(11)*y(53)/y(28)*y(11)/y(28);
  g1(17,41)=(-(1/y(41)));
  g1(17,42)=(-(y(1)*y(28)/y(8)*1/y(42)));
  g1(17,59)=y(50)*T191*1/y(59);
  g1(18,11)=(-(y(50)*params(2)*(params(7)-1)*y(48)/(params(7)*y(33)*y(26))));
  g1(18,48)=(-(y(50)*params(2)*y(11)*(params(7)-1)/(params(7)*y(33)*y(26))));
  g1(18,50)=(-(params(2)*y(11)*(params(7)-1)*y(48)/(params(7)*y(33)*y(26))+(1-params(5))*y(52)));
  g1(18,26)=(-(y(50)*(-(params(2)*y(11)*(params(7)-1)*y(48)*params(7)*y(33)))/(params(7)*y(33)*y(26)*params(7)*y(33)*y(26))));
  g1(18,27)=1;
  g1(18,52)=(-(y(50)*(1-params(5))));
  g1(18,33)=(-(y(50)*(-(params(2)*y(11)*(params(7)-1)*y(48)*params(7)*y(26)))/(params(7)*y(33)*y(26)*params(7)*y(33)*y(26))));
  g1(19,9)=(-(params(12)*1/y(9)));
  g1(19,29)=1/y(29);
  g1(19,60)=(-params(14));
  g1(19,10)=params(13)*1/y(10);
  g1(20,12)=(-y(14));
  g1(20,13)=(-y(14));
  g1(20,14)=(-(y(12)+y(13)-1));
  g1(20,15)=(-1);
  g1(20,26)=(-y(27));
  g1(20,27)=(-y(26));
  g1(20,30)=1;
  g1(20,31)=(-1);
  g1(21,11)=(-(y(50)*(y(45)*y(44)+y(54))));
  g1(21,44)=(-(y(45)*y(11)*y(50)));
  g1(21,45)=(-(y(11)*y(50)*y(44)));
  g1(21,50)=(-(y(11)*(y(45)*y(44)+y(54))));
  g1(21,31)=1;
  g1(21,54)=(-(y(11)*y(50)));
  g1(22,17)=(-1);
  g1(22,32)=1;
  g1(23,1)=(-(T249*T244*y(17)/y(4)*T358));
  g1(23,4)=(-(T249*T244*T358*T436));
  g1(23,17)=(-(T249*T244*T358*y(1)/y(4)));
  g1(23,39)=1;
  g1(24,1)=(-(T252*params(17)*y(17)/y(4)));
  g1(24,4)=(-(T252*params(17)*T436));
  g1(24,17)=(-(T252*params(17)*y(1)/y(4)));
  g1(24,40)=1;
  g1(25,1)=(-(T259*T255*y(28)/y(8)*T366));
  g1(25,8)=(-(T259*T255*T366*T605));
  g1(25,28)=(-(T259*T255*T366*y(1)/y(8)));
  g1(25,41)=1;
  g1(26,1)=(-(T262*params(18)*y(28)/y(8)));
  g1(26,8)=(-(T262*params(18)*T605));
  g1(26,28)=(-(T262*params(18)*y(1)/y(8)));
  g1(26,42)=1;
  g1(27,34)=getPowerDeriv(y(34),1-params(23),1);
  g1(27,35)=(-((1-params(24))*getPowerDeriv(y(35),1-params(23),1)));
  g1(28,34)=(-T280);
  g1(28,35)=1;
  g1(28,36)=(-(y(34)*params(23)/(params(23)-1)*1/y(37)));
  g1(28,37)=(-(y(34)*params(23)/(params(23)-1)*(-y(36))/(y(37)*y(37))));
  g1(29,11)=(-(y(56)*y(55)^params(23)*params(1)*params(24)*getPowerDeriv(y(11),1-params(4),1)));
  g1(29,19)=(-(T58*y(22)));
  g1(29,22)=(-(T58*y(19)));
  g1(29,33)=(-(y(19)*y(22)*T648));
  g1(29,55)=(-(y(56)*T288*getPowerDeriv(y(55),params(23),1)));
  g1(29,36)=1;
  g1(29,56)=(-(T288*y(55)^params(23)));
  g1(30,11)=(-(y(57)*y(55)^(params(23)-1)*params(1)*params(24)*getPowerDeriv(y(11),1-params(4),1)));
  g1(30,19)=(-y(22));
  g1(30,22)=(-y(19));
  g1(30,55)=(-(y(57)*T288*getPowerDeriv(y(55),params(23)-1,1)));
  g1(30,37)=1;
  g1(30,57)=(-(T288*y(55)^(params(23)-1)));
  g1(31,50)=(-(y(38)/y(55)));
  g1(31,55)=(-((-(y(50)*y(38)))/(y(55)*y(55))));
  g1(31,38)=(-(y(50)/y(55)));
  g1(32,33)=(-(y(34)^params(21)*T648/(1/params(16))*getPowerDeriv(T313,params(22),1)));
  g1(32,34)=(-(T313^params(22)*getPowerDeriv(y(34),params(21),1)));
  g1(32,38)=1/params(20);
  g1(33,9)=(-1);
  g1(33,43)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],33,3600);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],33,216000);
end
end
