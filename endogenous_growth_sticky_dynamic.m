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

residual = zeros(35, 1);
T34 = (y(1)*y(19)/y(7))^params(10);
T60 = 1/y(35);
T69 = y(7)/y(1);
T71 = 1/T69^params(10);
T74 = y(19)^(1-params(10));
T75 = 1/T74;
T92 = (y(11)*y(49)/y(19))^2;
T115 = y(26)^(1-params(2));
T139 = y(1)^(-params(4));
T146 = params(6)/(1+params(3));
T148 = y(26)^(1+params(3));
T150 = y(22)-y(27)*T146*T148;
T151 = T150^(-params(4));
T160 = (y(6)/(y(1)*y(22)))^(1-params(8));
T167 = y(11)^(-params(4))*y(53);
T171 = (y(11)*y(51)/y(27))^params(8);
T182 = y(27)*params(6)*y(26)^params(3)*1/y(24);
T187 = (1-params(2))*T60*(params(7)-1)/params(7);
T193 = (y(6)/y(1))^(1-params(8));
T213 = (y(11)*y(55)/y(30))^2;
T268 = params(16)/2;
T273 = exp(T268*(y(1)*y(19)/y(4)-params(18))^2);
T276 = exp(params(16)*(y(1)*y(19)/y(4)-params(18)));
T279 = params(17)/2;
T283 = exp(T279*(y(1)*y(30)/y(8)-params(18))^2);
T286 = exp(params(17)*(y(1)*y(30)/y(8)-params(18)));
T304 = params(22)/(params(22)-1)*y(38)/y(39);
T312 = params(1)*params(23)*y(11)^(1-params(4));
T337 = T60/(1/params(15));
T351 = getPowerDeriv(y(1)*y(19)/y(7),params(10),1);
T361 = T69^params(10)*T69^params(10);
T379 = getPowerDeriv(y(6)/(y(1)*y(22)),1-params(8),1);
T385 = getPowerDeriv(y(6)/y(1),1-params(8),1);
T394 = 2*(y(1)*y(19)/y(4)-params(18));
T402 = 2*(y(1)*y(30)/y(8)-params(18));
T418 = getPowerDeriv(y(11)*y(51)/y(27),params(8),1);
T482 = (-(y(1)*y(19)))/(y(4)*y(4));
T532 = getPowerDeriv(T150,(-params(4)),1);
T599 = T532*(-(y(27)*T146*getPowerDeriv(y(26),1+params(3),1)));
T620 = 1/y(1);
T671 = (-(y(1)*y(30)))/(y(8)*y(8));
T716 = (-1)/(y(35)*y(35));
lhs =y(11);
rhs =params(9)+y(13)*params(9)*(y(12)-1);
residual(1)= lhs-rhs;
lhs =y(12)*y(1);
rhs =params(9)*y(2)+y(3);
residual(2)= lhs-rhs;
lhs =y(14);
rhs =y(12)*params(14)*y(31)*T34;
residual(3)= lhs-rhs;
lhs =y(16);
rhs =(-y(15))+params(9)*y(52)*(y(13)*y(48)+(1-y(13))*y(47));
residual(4)= lhs-rhs;
lhs =y(17);
rhs =y(18)+params(9)*y(52)*y(48);
residual(5)= lhs-rhs;
lhs =y(18);
rhs =1/params(7)*T60*y(20);
residual(6)= lhs-rhs;
lhs =y(12)*y(31)*params(14)*y(52)*y(47)*T71*T75;
rhs =1+log(y(42))*y(1)*y(19)/y(4)+log(y(41))-y(52)*log(y(60))*T92;
residual(7)= lhs-rhs;
lhs =y(52)*params(9)*params(25)*params(24)*(y(48)-y(47));
rhs =y(15)^(1-params(25));
residual(8)= lhs-rhs;
lhs =y(13);
rhs =params(24)*y(15)^params(25);
residual(9)= lhs-rhs;
lhs =y(21);
rhs =y(20);
residual(10)= lhs-rhs;
lhs =y(20);
rhs =T69^params(2)*T115;
residual(11)= lhs-rhs;
lhs =y(21);
rhs =y(22)+(1+log(y(43)))*y(30)+y(19)*(1+log(y(41)))+(y(12)-1)*y(15);
residual(12)= lhs-rhs;
lhs =y(23);
rhs =params(1)*y(24)/y(5)*T139;
residual(13)= lhs-rhs;
lhs =y(24);
rhs =T151+(-y(25))*params(8)*T160;
residual(14)= lhs-rhs;
lhs =y(25);
rhs =params(1)*(1-params(8))*T167*T171+T148*T146*T151;
residual(15)= lhs-rhs;
lhs =T151*T182;
rhs =T187*y(21)/y(26);
residual(16)= lhs-rhs;
lhs =y(27);
rhs =y(22)^params(8)*T193;
residual(17)= lhs-rhs;
lhs =y(28);
rhs =y(30)+T69*(1-params(5));
residual(18)= lhs-rhs;
lhs =y(29);
rhs =1+log(y(43))+y(1)*y(30)/y(8)*log(y(44))-y(52)*T213*log(y(61));
residual(19)= lhs-rhs;
lhs =y(29);
rhs =y(52)*(params(2)*y(11)*(params(7)-1)*y(50)/(params(7)*y(35)*y(28))+(1-params(5))*y(54));
residual(20)= lhs-rhs;
lhs =log(y(31));
rhs =params(13)*x(it_, 1)+params(11)*log(y(9))-params(12)*log(y(10));
residual(21)= lhs-rhs;
lhs =y(32);
rhs =y(17)+y(28)*y(29)+y(16)*(y(12)+y(14)-1)+y(33);
residual(22)= lhs-rhs;
lhs =y(33);
rhs =y(11)*y(52)*(y(47)*y(46)+y(56));
residual(23)= lhs-rhs;
lhs =y(34);
rhs =y(19);
residual(24)= lhs-rhs;
lhs =y(41);
rhs =T273;
residual(25)= lhs-rhs;
lhs =y(42);
rhs =T276;
residual(26)= lhs-rhs;
lhs =y(43);
rhs =T283;
residual(27)= lhs-rhs;
lhs =y(44);
rhs =T286;
residual(28)= lhs-rhs;
lhs =y(36)^(1-params(22));
rhs =params(23)+(1-params(23))*y(37)^(1-params(22));
residual(29)= lhs-rhs;
lhs =y(37);
rhs =y(36)*T304;
residual(30)= lhs-rhs;
lhs =y(38);
rhs =y(21)*T60*y(24)+T312*y(57)^params(22)*y(58);
residual(31)= lhs-rhs;
lhs =y(39);
rhs =y(21)*y(24)+T312*y(57)^(params(22)-1)*y(59);
residual(32)= lhs-rhs;
lhs =1;
rhs =y(52)*y(40)/y(57);
residual(33)= lhs-rhs;
lhs =y(40)/params(19);
rhs =y(36)^params(20)*T337^params(21);
residual(34)= lhs-rhs;
lhs =y(45);
rhs =y(9);
residual(35)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(35, 62);

  %
  % Jacobian matrix
  %

  g1(1,11)=1;
  g1(1,12)=(-(y(13)*params(9)));
  g1(1,13)=(-(params(9)*(y(12)-1)));
  g1(2,1)=y(12);
  g1(2,2)=(-params(9));
  g1(2,12)=y(1);
  g1(2,3)=(-1);
  g1(3,1)=(-(y(12)*params(14)*y(31)*y(19)/y(7)*T351));
  g1(3,12)=(-(params(14)*y(31)*T34));
  g1(3,14)=1;
  g1(3,19)=(-(y(12)*params(14)*y(31)*T351*y(1)/y(7)));
  g1(3,7)=(-(y(12)*params(14)*y(31)*T351*(-(y(1)*y(19)))/(y(7)*y(7))));
  g1(3,31)=(-(T34*y(12)*params(14)));
  g1(4,13)=(-(params(9)*y(52)*(y(48)-y(47))));
  g1(4,15)=1;
  g1(4,16)=1;
  g1(4,47)=(-(params(9)*y(52)*(1-y(13))));
  g1(4,48)=(-(y(13)*params(9)*y(52)));
  g1(4,52)=(-(params(9)*(y(13)*y(48)+(1-y(13))*y(47))));
  g1(5,17)=1;
  g1(5,48)=(-(params(9)*y(52)));
  g1(5,18)=(-1);
  g1(5,52)=(-(params(9)*y(48)));
  g1(6,18)=1;
  g1(6,20)=(-(1/params(7)*T60));
  g1(6,35)=(-(y(20)*1/params(7)*T716));
  g1(7,1)=T75*y(12)*y(31)*params(14)*y(52)*y(47)*(-((-y(7))/(y(1)*y(1))*getPowerDeriv(T69,params(10),1)))/T361-log(y(42))*y(19)/y(4);
  g1(7,11)=y(52)*log(y(60))*y(49)/y(19)*2*y(11)*y(49)/y(19);
  g1(7,12)=T75*y(31)*params(14)*y(52)*y(47)*T71;
  g1(7,47)=T75*T71*y(12)*y(31)*params(14)*y(52);
  g1(7,4)=(-(log(y(42))*T482));
  g1(7,19)=y(12)*y(31)*params(14)*y(52)*y(47)*T71*(-(getPowerDeriv(y(19),1-params(10),1)))/(T74*T74)-(log(y(42))*y(1)/y(4)-y(52)*log(y(60))*2*y(11)*y(49)/y(19)*(-(y(11)*y(49)))/(y(19)*y(19)));
  g1(7,49)=y(52)*log(y(60))*2*y(11)*y(49)/y(19)*y(11)/y(19);
  g1(7,52)=T75*T71*y(12)*y(31)*params(14)*y(47)-(-(log(y(60))*T92));
  g1(7,7)=T75*y(12)*y(31)*params(14)*y(52)*y(47)*(-(getPowerDeriv(T69,params(10),1)*T620))/T361;
  g1(7,31)=T75*T71*y(12)*params(14)*y(52)*y(47);
  g1(7,41)=(-(1/y(41)));
  g1(7,42)=(-(y(1)*y(19)/y(4)*1/y(42)));
  g1(7,60)=T92*y(52)*1/y(60);
  g1(8,15)=(-(getPowerDeriv(y(15),1-params(25),1)));
  g1(8,47)=(-(y(52)*params(9)*params(25)*params(24)));
  g1(8,48)=y(52)*params(9)*params(25)*params(24);
  g1(8,52)=params(9)*params(25)*params(24)*(y(48)-y(47));
  g1(9,13)=1;
  g1(9,15)=(-(params(24)*getPowerDeriv(y(15),params(25),1)));
  g1(10,20)=(-1);
  g1(10,21)=1;
  g1(11,1)=(-(T115*(-y(7))/(y(1)*y(1))*getPowerDeriv(T69,params(2),1)));
  g1(11,20)=1;
  g1(11,26)=(-(T69^params(2)*getPowerDeriv(y(26),1-params(2),1)));
  g1(11,7)=(-(T115*getPowerDeriv(T69,params(2),1)*T620));
  g1(12,12)=(-y(15));
  g1(12,15)=(-(y(12)-1));
  g1(12,19)=(-(1+log(y(41))));
  g1(12,21)=1;
  g1(12,22)=(-1);
  g1(12,30)=(-(1+log(y(43))));
  g1(12,41)=(-(y(19)*1/y(41)));
  g1(12,43)=(-(y(30)*1/y(43)));
  g1(13,1)=(-(params(1)*y(24)/y(5)*getPowerDeriv(y(1),(-params(4)),1)));
  g1(13,23)=1;
  g1(13,5)=(-(T139*(-(params(1)*y(24)))/(y(5)*y(5))));
  g1(13,24)=(-(T139*params(1)/y(5)));
  g1(14,1)=(-((-y(25))*params(8)*(-(y(22)*y(6)))/(y(1)*y(22)*y(1)*y(22))*T379));
  g1(14,22)=(-(T532+(-y(25))*params(8)*T379*(-(y(1)*y(6)))/(y(1)*y(22)*y(1)*y(22))));
  g1(14,24)=1;
  g1(14,25)=(-(T160*(-params(8))));
  g1(14,26)=(-T599);
  g1(14,6)=(-((-y(25))*params(8)*T379*1/(y(1)*y(22))));
  g1(14,27)=(-(T532*(-(T146*T148))));
  g1(15,11)=(-(params(1)*(1-params(8))*(T171*y(53)*getPowerDeriv(y(11),(-params(4)),1)+T167*y(51)/y(27)*T418)));
  g1(15,22)=(-(T148*T146*T532));
  g1(15,51)=(-(params(1)*(1-params(8))*T167*T418*y(11)/y(27)));
  g1(15,25)=1;
  g1(15,53)=(-(params(1)*(1-params(8))*y(11)^(-params(4))*T171));
  g1(15,26)=(-(T146*T151*getPowerDeriv(y(26),1+params(3),1)+T148*T146*T599));
  g1(15,27)=(-(params(1)*(1-params(8))*T167*T418*(-(y(11)*y(51)))/(y(27)*y(27))+T148*T146*T532*(-(T146*T148))));
  g1(16,21)=(-(T187*1/y(26)));
  g1(16,22)=T182*T532;
  g1(16,24)=T151*y(27)*params(6)*y(26)^params(3)*(-1)/(y(24)*y(24));
  g1(16,26)=T182*T599+T151*1/y(24)*y(27)*params(6)*getPowerDeriv(y(26),params(3),1)-T187*(-y(21))/(y(26)*y(26));
  g1(16,27)=T182*T532*(-(T146*T148))+T151*1/y(24)*params(6)*y(26)^params(3);
  g1(16,35)=(-(y(21)/y(26)*(1-params(2))*(params(7)-1)/params(7)*T716));
  g1(17,1)=(-(y(22)^params(8)*(-y(6))/(y(1)*y(1))*T385));
  g1(17,22)=(-(T193*getPowerDeriv(y(22),params(8),1)));
  g1(17,6)=(-(y(22)^params(8)*T385*T620));
  g1(17,27)=1;
  g1(18,1)=(-((1-params(5))*(-y(7))/(y(1)*y(1))));
  g1(18,7)=(-((1-params(5))*T620));
  g1(18,28)=1;
  g1(18,30)=(-1);
  g1(19,1)=(-(log(y(44))*y(30)/y(8)));
  g1(19,11)=log(y(61))*y(52)*y(55)/y(30)*2*y(11)*y(55)/y(30);
  g1(19,52)=T213*log(y(61));
  g1(19,29)=1;
  g1(19,8)=(-(log(y(44))*T671));
  g1(19,30)=(-(log(y(44))*y(1)/y(8)-log(y(61))*y(52)*2*y(11)*y(55)/y(30)*(-(y(11)*y(55)))/(y(30)*y(30))));
  g1(19,55)=log(y(61))*y(52)*2*y(11)*y(55)/y(30)*y(11)/y(30);
  g1(19,43)=(-(1/y(43)));
  g1(19,44)=(-(y(1)*y(30)/y(8)*1/y(44)));
  g1(19,61)=y(52)*T213*1/y(61);
  g1(20,11)=(-(y(52)*params(2)*(params(7)-1)*y(50)/(params(7)*y(35)*y(28))));
  g1(20,50)=(-(y(52)*params(2)*y(11)*(params(7)-1)/(params(7)*y(35)*y(28))));
  g1(20,52)=(-(params(2)*y(11)*(params(7)-1)*y(50)/(params(7)*y(35)*y(28))+(1-params(5))*y(54)));
  g1(20,28)=(-(y(52)*(-(params(2)*y(11)*(params(7)-1)*y(50)*params(7)*y(35)))/(params(7)*y(35)*y(28)*params(7)*y(35)*y(28))));
  g1(20,29)=1;
  g1(20,54)=(-(y(52)*(1-params(5))));
  g1(20,35)=(-(y(52)*(-(params(2)*y(11)*(params(7)-1)*y(50)*params(7)*y(28)))/(params(7)*y(35)*y(28)*params(7)*y(35)*y(28))));
  g1(21,9)=(-(params(11)*1/y(9)));
  g1(21,31)=1/y(31);
  g1(21,62)=(-params(13));
  g1(21,10)=params(12)*1/y(10);
  g1(22,12)=(-y(16));
  g1(22,14)=(-y(16));
  g1(22,16)=(-(y(12)+y(14)-1));
  g1(22,17)=(-1);
  g1(22,28)=(-y(29));
  g1(22,29)=(-y(28));
  g1(22,32)=1;
  g1(22,33)=(-1);
  g1(23,11)=(-(y(52)*(y(47)*y(46)+y(56))));
  g1(23,46)=(-(y(47)*y(11)*y(52)));
  g1(23,47)=(-(y(11)*y(52)*y(46)));
  g1(23,52)=(-(y(11)*(y(47)*y(46)+y(56))));
  g1(23,33)=1;
  g1(23,56)=(-(y(11)*y(52)));
  g1(24,19)=(-1);
  g1(24,34)=1;
  g1(25,1)=(-(T273*T268*y(19)/y(4)*T394));
  g1(25,4)=(-(T273*T268*T394*T482));
  g1(25,19)=(-(T273*T268*T394*y(1)/y(4)));
  g1(25,41)=1;
  g1(26,1)=(-(T276*params(16)*y(19)/y(4)));
  g1(26,4)=(-(T276*params(16)*T482));
  g1(26,19)=(-(T276*params(16)*y(1)/y(4)));
  g1(26,42)=1;
  g1(27,1)=(-(T283*T279*y(30)/y(8)*T402));
  g1(27,8)=(-(T283*T279*T402*T671));
  g1(27,30)=(-(T283*T279*T402*y(1)/y(8)));
  g1(27,43)=1;
  g1(28,1)=(-(T286*params(17)*y(30)/y(8)));
  g1(28,8)=(-(T286*params(17)*T671));
  g1(28,30)=(-(T286*params(17)*y(1)/y(8)));
  g1(28,44)=1;
  g1(29,36)=getPowerDeriv(y(36),1-params(22),1);
  g1(29,37)=(-((1-params(23))*getPowerDeriv(y(37),1-params(22),1)));
  g1(30,36)=(-T304);
  g1(30,37)=1;
  g1(30,38)=(-(y(36)*params(22)/(params(22)-1)*1/y(39)));
  g1(30,39)=(-(y(36)*params(22)/(params(22)-1)*(-y(38))/(y(39)*y(39))));
  g1(31,11)=(-(y(58)*y(57)^params(22)*params(1)*params(23)*getPowerDeriv(y(11),1-params(4),1)));
  g1(31,21)=(-(T60*y(24)));
  g1(31,24)=(-(T60*y(21)));
  g1(31,35)=(-(y(21)*y(24)*T716));
  g1(31,57)=(-(y(58)*T312*getPowerDeriv(y(57),params(22),1)));
  g1(31,38)=1;
  g1(31,58)=(-(T312*y(57)^params(22)));
  g1(32,11)=(-(y(59)*y(57)^(params(22)-1)*params(1)*params(23)*getPowerDeriv(y(11),1-params(4),1)));
  g1(32,21)=(-y(24));
  g1(32,24)=(-y(21));
  g1(32,57)=(-(y(59)*T312*getPowerDeriv(y(57),params(22)-1,1)));
  g1(32,39)=1;
  g1(32,59)=(-(T312*y(57)^(params(22)-1)));
  g1(33,52)=(-(y(40)/y(57)));
  g1(33,57)=(-((-(y(52)*y(40)))/(y(57)*y(57))));
  g1(33,40)=(-(y(52)/y(57)));
  g1(34,35)=(-(y(36)^params(20)*T716/(1/params(15))*getPowerDeriv(T337,params(21),1)));
  g1(34,36)=(-(T337^params(21)*getPowerDeriv(y(36),params(20),1)));
  g1(34,40)=1/params(19);
  g1(35,9)=(-1);
  g1(35,45)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],35,3844);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],35,238328);
end
end
