function [residual, g1, g2] = endogenous_growth_sticky_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                    columns: equations in order of declaration
%                                                    rows: variables in declaration order
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: equations in order of declaration
%                                                       rows: variables in declaration order
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 35, 1);

%
% Model equations
%

T31 = (y(1)*y(9)/y(18))^params(10);
T55 = 1/y(25);
T64 = y(18)/y(1);
T66 = 1/T64^params(10);
T69 = y(9)^(1-params(10));
T70 = 1/T69;
T74 = y(1)*y(9)/y(9);
T104 = y(16)^(1-params(2));
T126 = y(1)^(-params(4));
T133 = params(6)/(1+params(3));
T135 = y(16)^(1+params(3));
T137 = y(12)-y(17)*T133*T135;
T138 = T137^(-params(4));
T146 = (y(17)/(y(1)*y(12)))^(1-params(8));
T151 = T126*y(15);
T153 = (y(1)*y(12)/y(17))^params(8);
T164 = y(17)*params(6)*y(16)^params(3)*1/y(14);
T169 = (1-params(2))*T55*(params(7)-1)/params(7);
T175 = (y(17)/y(1))^(1-params(8));
T185 = y(1)*y(20)/y(20);
T241 = exp(params(16)/2*(T74-params(18))^2);
T251 = exp(params(17)/2*(T185-params(18))^2);
T272 = params(22)/(params(22)-1)*y(28)/y(29);
T280 = params(1)*params(23)*y(1)^(1-params(4));
T302 = T55/(1/params(15));
T309 = getPowerDeriv(y(1)*y(9)/y(18),params(10),1);
T319 = T64^params(10)*T64^params(10);
T338 = getPowerDeriv(y(17)/(y(1)*y(12)),1-params(8),1);
T344 = getPowerDeriv(y(1)*y(12)/y(17),params(8),1);
T353 = getPowerDeriv(y(17)/y(1),1-params(8),1);
T445 = getPowerDeriv(T137,(-params(4)),1);
T504 = T445*(-(y(17)*T133*getPowerDeriv(y(16),1+params(3),1)));
T544 = 1/y(1);
T587 = (-1)/(y(25)*y(25));
lhs =y(1);
rhs =params(9)+params(9)*y(3)*(y(2)-1);
residual(1)= lhs-rhs;
lhs =y(1)*y(2);
rhs =params(9)*y(2)+y(4);
residual(2)= lhs-rhs;
lhs =y(4);
rhs =y(2)*params(14)*y(21)*T31;
residual(3)= lhs-rhs;
lhs =y(6);
rhs =(-y(5))+params(9)*y(13)*(y(3)*y(7)+y(6)*(1-y(3)));
residual(4)= lhs-rhs;
lhs =y(7);
rhs =y(8)+params(9)*y(13)*y(7);
residual(5)= lhs-rhs;
lhs =y(8);
rhs =1/params(7)*T55*y(10);
residual(6)= lhs-rhs;
lhs =y(2)*y(21)*params(14)*y(6)*y(13)*T66*T70;
rhs =1+log(y(32))*T74+log(y(31))-y(13)*log(y(32))*T74^2;
residual(7)= lhs-rhs;
lhs =y(13)*params(9)*params(25)*params(24)*(y(7)-y(6));
rhs =y(5)^(1-params(25));
residual(8)= lhs-rhs;
lhs =y(3);
rhs =params(24)*y(5)^params(25);
residual(9)= lhs-rhs;
lhs =y(11);
rhs =y(10);
residual(10)= lhs-rhs;
lhs =y(10);
rhs =T64^params(2)*T104;
residual(11)= lhs-rhs;
lhs =y(11);
rhs =y(12)+(1+log(y(33)))*y(20)+y(9)*(1+log(y(31)))+(y(2)-1)*y(5);
residual(12)= lhs-rhs;
lhs =y(13);
rhs =params(1)*y(14)/y(14)*T126;
residual(13)= lhs-rhs;
lhs =y(14);
rhs =T138+(-y(15))*params(8)*T146;
residual(14)= lhs-rhs;
lhs =y(15);
rhs =params(1)*(1-params(8))*T151*T153+T135*T133*T138;
residual(15)= lhs-rhs;
lhs =T138*T164;
rhs =T169*y(11)/y(16);
residual(16)= lhs-rhs;
lhs =y(17);
rhs =y(12)^params(8)*T175;
residual(17)= lhs-rhs;
lhs =y(18);
rhs =y(20)+T64*(1-params(5));
residual(18)= lhs-rhs;
lhs =y(19);
rhs =1+log(y(33))+T185*log(y(34))-log(y(34))*y(13)*T185^2;
residual(19)= lhs-rhs;
lhs =y(19);
rhs =y(13)*(params(2)*y(10)*y(1)*(params(7)-1)/(params(7)*y(18)*y(25))+(1-params(5))*y(19));
residual(20)= lhs-rhs;
lhs =log(y(21));
rhs =params(13)*x(1)+log(y(21))*params(11)-params(12)*log(y(35));
residual(21)= lhs-rhs;
lhs =y(22);
rhs =y(7)+y(18)*y(19)+y(6)*(y(2)+y(4)-1)+y(23);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =y(1)*y(13)*(y(23)+y(4)*y(6));
residual(23)= lhs-rhs;
lhs =y(24);
rhs =y(9);
residual(24)= lhs-rhs;
lhs =y(31);
rhs =T241;
residual(25)= lhs-rhs;
lhs =y(32);
rhs =exp(params(16)*(T74-params(18)));
residual(26)= lhs-rhs;
lhs =y(33);
rhs =T251;
residual(27)= lhs-rhs;
lhs =y(34);
rhs =exp(params(17)*(T185-params(18)));
residual(28)= lhs-rhs;
lhs =y(26)^(1-params(22));
rhs =params(23)+(1-params(23))*y(27)^(1-params(22));
residual(29)= lhs-rhs;
lhs =y(27);
rhs =y(26)*T272;
residual(30)= lhs-rhs;
lhs =y(28);
rhs =y(11)*T55*y(14)+y(28)*T280*y(26)^params(22);
residual(31)= lhs-rhs;
lhs =y(29);
rhs =y(11)*y(14)+y(29)*T280*y(26)^(params(22)-1);
residual(32)= lhs-rhs;
lhs =1;
rhs =y(13)*y(30)/y(26);
residual(33)= lhs-rhs;
lhs =y(30)/params(19);
rhs =y(26)^params(20)*T302^params(21);
residual(34)= lhs-rhs;
lhs =y(35);
rhs =y(21);
residual(35)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(35, 35);

  %
  % Jacobian matrix
  %

  g1(1,1)=1;
  g1(1,2)=(-(params(9)*y(3)));
  g1(1,3)=(-(params(9)*(y(2)-1)));
  g1(2,1)=y(2);
  g1(2,2)=y(1)-params(9);
  g1(2,4)=(-1);
  g1(3,1)=(-(y(2)*params(14)*y(21)*y(9)/y(18)*T309));
  g1(3,2)=(-(params(14)*y(21)*T31));
  g1(3,4)=1;
  g1(3,9)=(-(y(2)*params(14)*y(21)*T309*y(1)/y(18)));
  g1(3,18)=(-(y(2)*params(14)*y(21)*T309*(-(y(1)*y(9)))/(y(18)*y(18))));
  g1(3,21)=(-(T31*y(2)*params(14)));
  g1(4,3)=(-(params(9)*y(13)*(y(7)-y(6))));
  g1(4,5)=1;
  g1(4,6)=1-params(9)*y(13)*(1-y(3));
  g1(4,7)=(-(y(3)*params(9)*y(13)));
  g1(4,13)=(-(params(9)*(y(3)*y(7)+y(6)*(1-y(3)))));
  g1(5,7)=1-params(9)*y(13);
  g1(5,8)=(-1);
  g1(5,13)=(-(params(9)*y(7)));
  g1(6,8)=1;
  g1(6,10)=(-(1/params(7)*T55));
  g1(6,25)=(-(y(10)*1/params(7)*T587));
  g1(7,1)=T70*y(2)*y(21)*params(14)*y(6)*y(13)*(-((-y(18))/(y(1)*y(1))*getPowerDeriv(T64,params(10),1)))/T319-(log(y(32))-y(13)*log(y(32))*2*T74);
  g1(7,2)=T70*y(21)*params(14)*y(6)*y(13)*T66;
  g1(7,6)=T70*T66*y(2)*y(21)*params(14)*y(13);
  g1(7,9)=y(2)*y(21)*params(14)*y(6)*y(13)*T66*(-(getPowerDeriv(y(9),1-params(10),1)))/(T69*T69);
  g1(7,13)=T70*T66*y(2)*y(21)*params(14)*y(6)-(-(log(y(32))*T74^2));
  g1(7,18)=T70*y(2)*y(21)*params(14)*y(6)*y(13)*(-(getPowerDeriv(T64,params(10),1)*T544))/T319;
  g1(7,21)=T70*T66*y(2)*params(14)*y(6)*y(13);
  g1(7,31)=(-(1/y(31)));
  g1(7,32)=(-(T74*1/y(32)-T74^2*y(13)*1/y(32)));
  g1(8,5)=(-(getPowerDeriv(y(5),1-params(25),1)));
  g1(8,6)=(-(y(13)*params(9)*params(25)*params(24)));
  g1(8,7)=y(13)*params(9)*params(25)*params(24);
  g1(8,13)=params(9)*params(25)*params(24)*(y(7)-y(6));
  g1(9,3)=1;
  g1(9,5)=(-(params(24)*getPowerDeriv(y(5),params(25),1)));
  g1(10,10)=(-1);
  g1(10,11)=1;
  g1(11,1)=(-(T104*(-y(18))/(y(1)*y(1))*getPowerDeriv(T64,params(2),1)));
  g1(11,10)=1;
  g1(11,16)=(-(T64^params(2)*getPowerDeriv(y(16),1-params(2),1)));
  g1(11,18)=(-(T104*getPowerDeriv(T64,params(2),1)*T544));
  g1(12,2)=(-y(5));
  g1(12,5)=(-(y(2)-1));
  g1(12,9)=(-(1+log(y(31))));
  g1(12,11)=1;
  g1(12,12)=(-1);
  g1(12,20)=(-(1+log(y(33))));
  g1(12,31)=(-(y(9)*1/y(31)));
  g1(12,33)=(-(y(20)*1/y(33)));
  g1(13,1)=(-(params(1)*y(14)/y(14)*getPowerDeriv(y(1),(-params(4)),1)));
  g1(13,13)=1;
  g1(14,1)=(-((-y(15))*params(8)*(-(y(12)*y(17)))/(y(1)*y(12)*y(1)*y(12))*T338));
  g1(14,12)=(-(T445+(-y(15))*params(8)*T338*(-(y(1)*y(17)))/(y(1)*y(12)*y(1)*y(12))));
  g1(14,14)=1;
  g1(14,15)=(-(T146*(-params(8))));
  g1(14,16)=(-T504);
  g1(14,17)=(-(T445*(-(T133*T135))+(-y(15))*params(8)*T338*1/(y(1)*y(12))));
  g1(15,1)=(-(params(1)*(1-params(8))*(T153*y(15)*getPowerDeriv(y(1),(-params(4)),1)+T151*y(12)/y(17)*T344)));
  g1(15,12)=(-(params(1)*(1-params(8))*T151*T344*y(1)/y(17)+T135*T133*T445));
  g1(15,15)=1-params(1)*(1-params(8))*T126*T153;
  g1(15,16)=(-(T133*T138*getPowerDeriv(y(16),1+params(3),1)+T135*T133*T504));
  g1(15,17)=(-(params(1)*(1-params(8))*T151*T344*(-(y(1)*y(12)))/(y(17)*y(17))+T135*T133*T445*(-(T133*T135))));
  g1(16,11)=(-(T169*1/y(16)));
  g1(16,12)=T164*T445;
  g1(16,14)=T138*y(17)*params(6)*y(16)^params(3)*(-1)/(y(14)*y(14));
  g1(16,16)=T164*T504+T138*1/y(14)*y(17)*params(6)*getPowerDeriv(y(16),params(3),1)-T169*(-y(11))/(y(16)*y(16));
  g1(16,17)=T164*T445*(-(T133*T135))+T138*1/y(14)*params(6)*y(16)^params(3);
  g1(16,25)=(-(y(11)/y(16)*(1-params(2))*(params(7)-1)/params(7)*T587));
  g1(17,1)=(-(y(12)^params(8)*(-y(17))/(y(1)*y(1))*T353));
  g1(17,12)=(-(T175*getPowerDeriv(y(12),params(8),1)));
  g1(17,17)=1-y(12)^params(8)*T353*T544;
  g1(18,1)=(-((1-params(5))*(-y(18))/(y(1)*y(1))));
  g1(18,18)=1-(1-params(5))*T544;
  g1(18,20)=(-1);
  g1(19,1)=(-(log(y(34))-log(y(34))*y(13)*2*T185));
  g1(19,13)=log(y(34))*T185^2;
  g1(19,19)=1;
  g1(19,33)=(-(1/y(33)));
  g1(19,34)=(-(T185*1/y(34)-y(13)*T185^2*1/y(34)));
  g1(20,1)=(-(y(13)*params(2)*y(10)*(params(7)-1)/(params(7)*y(18)*y(25))));
  g1(20,10)=(-(y(13)*params(2)*y(1)*(params(7)-1)/(params(7)*y(18)*y(25))));
  g1(20,13)=(-(params(2)*y(10)*y(1)*(params(7)-1)/(params(7)*y(18)*y(25))+(1-params(5))*y(19)));
  g1(20,18)=(-(y(13)*(-(params(2)*y(10)*y(1)*(params(7)-1)*params(7)*y(25)))/(params(7)*y(18)*y(25)*params(7)*y(18)*y(25))));
  g1(20,19)=1-y(13)*(1-params(5));
  g1(20,25)=(-(y(13)*(-(params(2)*y(10)*y(1)*(params(7)-1)*y(18)*params(7)))/(params(7)*y(18)*y(25)*params(7)*y(18)*y(25))));
  g1(21,21)=1/y(21)-params(11)*1/y(21);
  g1(21,35)=params(12)*1/y(35);
  g1(22,2)=(-y(6));
  g1(22,4)=(-y(6));
  g1(22,6)=(-(y(2)+y(4)-1));
  g1(22,7)=(-1);
  g1(22,18)=(-y(19));
  g1(22,19)=(-y(18));
  g1(22,22)=1;
  g1(22,23)=(-1);
  g1(23,1)=(-(y(13)*(y(23)+y(4)*y(6))));
  g1(23,4)=(-(y(6)*y(1)*y(13)));
  g1(23,6)=(-(y(4)*y(1)*y(13)));
  g1(23,13)=(-(y(1)*(y(23)+y(4)*y(6))));
  g1(23,23)=1-y(1)*y(13);
  g1(24,9)=(-1);
  g1(24,24)=1;
  g1(25,1)=(-(T241*params(16)/2*2*(T74-params(18))));
  g1(25,31)=1;
  g1(26,1)=(-(params(16)*exp(params(16)*(T74-params(18)))));
  g1(26,32)=1;
  g1(27,1)=(-(T251*params(17)/2*2*(T185-params(18))));
  g1(27,33)=1;
  g1(28,1)=(-(params(17)*exp(params(17)*(T185-params(18)))));
  g1(28,34)=1;
  g1(29,26)=getPowerDeriv(y(26),1-params(22),1);
  g1(29,27)=(-((1-params(23))*getPowerDeriv(y(27),1-params(22),1)));
  g1(30,26)=(-T272);
  g1(30,27)=1;
  g1(30,28)=(-(y(26)*params(22)/(params(22)-1)*1/y(29)));
  g1(30,29)=(-(y(26)*params(22)/(params(22)-1)*(-y(28))/(y(29)*y(29))));
  g1(31,1)=(-(y(28)*y(26)^params(22)*params(1)*params(23)*getPowerDeriv(y(1),1-params(4),1)));
  g1(31,11)=(-(T55*y(14)));
  g1(31,14)=(-(T55*y(11)));
  g1(31,25)=(-(y(11)*y(14)*T587));
  g1(31,26)=(-(y(28)*T280*getPowerDeriv(y(26),params(22),1)));
  g1(31,28)=1-T280*y(26)^params(22);
  g1(32,1)=(-(y(29)*y(26)^(params(22)-1)*params(1)*params(23)*getPowerDeriv(y(1),1-params(4),1)));
  g1(32,11)=(-y(14));
  g1(32,14)=(-y(11));
  g1(32,26)=(-(y(29)*T280*getPowerDeriv(y(26),params(22)-1,1)));
  g1(32,29)=1-T280*y(26)^(params(22)-1);
  g1(33,13)=(-(y(30)/y(26)));
  g1(33,26)=(-((-(y(13)*y(30)))/(y(26)*y(26))));
  g1(33,30)=(-(y(13)/y(26)));
  g1(34,25)=(-(y(26)^params(20)*T587/(1/params(15))*getPowerDeriv(T302,params(21),1)));
  g1(34,26)=(-(T302^params(21)*getPowerDeriv(y(26),params(20),1)));
  g1(34,30)=1/params(19);
  g1(35,21)=(-1);
  g1(35,35)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],35,1225);
end
end
