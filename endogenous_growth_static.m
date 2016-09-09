function [residual, g1, g2] = endogenous_growth_static(y, x, params)
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

residual = zeros( 26, 1);

%
% Model equations
%

T54 = 1/params(7)*1/params(15);
T62 = (y(2)/y(7))^(1-params(10));
T67 = y(1)*y(7)/y(7);
T83 = (y(16)/y(1))^params(2);
T86 = y(14)^(1-params(2));
T104 = y(1)^(-params(4));
T111 = params(6)/(1+params(3));
T113 = y(14)^(1+params(3));
T115 = y(10)-y(15)*T111*T113;
T116 = T115^(-params(4));
T124 = (y(15)/(y(1)*y(10)))^(1-params(8));
T129 = T104*y(13);
T131 = (y(1)*y(10)/y(15))^params(8);
T142 = y(15)*params(6)*y(14)^params(3)*1/y(12);
T147 = (1-params(2))*1/params(15)*(params(7)-1)/params(7);
T153 = (y(15)/y(1))^(1-params(8));
T163 = y(1)*y(18)/y(18);
T212 = exp(params(16)/2*(T67-params(18))^2);
T222 = exp(params(17)/2*(T163-params(18))^2);
T234 = getPowerDeriv(y(16)/y(1),params(2),1);
T245 = getPowerDeriv(y(15)/(y(1)*y(10)),1-params(8),1);
T251 = getPowerDeriv(y(1)*y(10)/y(15),params(8),1);
T260 = getPowerDeriv(y(15)/y(1),1-params(8),1);
T298 = getPowerDeriv(y(2)/y(7),1-params(10),1);
T330 = getPowerDeriv(T115,(-params(4)),1);
T377 = T330*(-(y(15)*T111*getPowerDeriv(y(14),1+params(3),1)));
lhs =y(1);
rhs =params(9)*(1+params(11)*(y(2)+y(3)-1));
residual(1)= lhs-rhs;
lhs =y(1)*y(2);
rhs =params(9)*(y(2)+y(3)*(1-params(11)));
residual(2)= lhs-rhs;
lhs =y(3);
rhs =params(14)*y(19)*y(2)^(1-params(10))*y(7)^params(10);
residual(3)= lhs-rhs;
lhs =y(4);
rhs =params(11)*y(5)+y(4)*params(9)*(1-params(11))*y(11);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =y(6)+params(9)*y(5)*y(11);
residual(5)= lhs-rhs;
lhs =y(6);
rhs =T54*y(8);
residual(6)= lhs-rhs;
lhs =y(19)*y(4)*params(14)*params(10)*T62;
rhs =1+log(y(24))*T67+log(y(23))-y(11)*log(y(24))*T67^2;
residual(7)= lhs-rhs;
lhs =y(9);
rhs =y(8);
residual(8)= lhs-rhs;
lhs =y(8);
rhs =T83*T86;
residual(9)= lhs-rhs;
lhs =y(9);
rhs =y(7)+y(10)+(1+log(y(25)))*y(18);
residual(10)= lhs-rhs;
lhs =y(11);
rhs =params(1)*y(12)/y(12)*T104;
residual(11)= lhs-rhs;
lhs =y(12);
rhs =T116+(-y(13))*params(8)*T124;
residual(12)= lhs-rhs;
lhs =y(13);
rhs =params(1)*(1-params(8))*T129*T131+T113*T111*T116;
residual(13)= lhs-rhs;
lhs =T116*T142;
rhs =T147*y(9)/y(14);
residual(14)= lhs-rhs;
lhs =y(15);
rhs =y(10)^params(8)*T153;
residual(15)= lhs-rhs;
lhs =y(16);
rhs =y(18)+y(16)/y(1)*(1-params(5));
residual(16)= lhs-rhs;
lhs =y(17);
rhs =1+log(y(25))+T163*log(y(26))-log(y(26))*y(11)*T163^2;
residual(17)= lhs-rhs;
lhs =y(17);
rhs =y(11)*(params(2)*y(8)*y(1)*(params(7)-1)/(params(7)*params(15)*y(16))+(1-params(5))*y(17));
residual(18)= lhs-rhs;
lhs =log(y(19));
rhs =log(y(19))*params(12)+params(13)*x(1);
residual(19)= lhs-rhs;
lhs =y(20);
rhs =y(5)+y(16)*y(17)+(y(2)+y(3)-1)*y(4)+y(21);
residual(20)= lhs-rhs;
lhs =y(21);
rhs =y(1)*y(11)*(y(21)+y(3)*y(4));
residual(21)= lhs-rhs;
lhs =y(22);
rhs =y(7);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =T212;
residual(23)= lhs-rhs;
lhs =y(24);
rhs =exp(params(16)*(T67-params(18)));
residual(24)= lhs-rhs;
lhs =y(25);
rhs =T222;
residual(25)= lhs-rhs;
lhs =y(26);
rhs =exp(params(17)*(T163-params(18)));
residual(26)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(26, 26);

  %
  % Jacobian matrix
  %

  g1(1,1)=1;
  g1(1,2)=(-(params(9)*params(11)));
  g1(1,3)=(-(params(9)*params(11)));
  g1(2,1)=y(2);
  g1(2,2)=y(1)-params(9);
  g1(2,3)=(-(params(9)*(1-params(11))));
  g1(3,2)=(-(y(7)^params(10)*params(14)*y(19)*getPowerDeriv(y(2),1-params(10),1)));
  g1(3,3)=1;
  g1(3,7)=(-(params(14)*y(19)*y(2)^(1-params(10))*getPowerDeriv(y(7),params(10),1)));
  g1(3,19)=(-(y(7)^params(10)*params(14)*y(2)^(1-params(10))));
  g1(4,4)=1-params(9)*(1-params(11))*y(11);
  g1(4,5)=(-params(11));
  g1(4,11)=(-(y(4)*params(9)*(1-params(11))));
  g1(5,5)=1-params(9)*y(11);
  g1(5,6)=(-1);
  g1(5,11)=(-(params(9)*y(5)));
  g1(6,6)=1;
  g1(6,8)=(-T54);
  g1(7,1)=(-(log(y(24))-y(11)*log(y(24))*2*T67));
  g1(7,2)=y(19)*y(4)*params(14)*params(10)*1/y(7)*T298;
  g1(7,4)=T62*y(19)*params(14)*params(10);
  g1(7,7)=y(19)*y(4)*params(14)*params(10)*T298*(-y(2))/(y(7)*y(7));
  g1(7,11)=log(y(24))*T67^2;
  g1(7,19)=y(4)*params(14)*params(10)*T62;
  g1(7,23)=(-(1/y(23)));
  g1(7,24)=(-(T67*1/y(24)-T67^2*y(11)*1/y(24)));
  g1(8,8)=(-1);
  g1(8,9)=1;
  g1(9,1)=(-(T86*(-y(16))/(y(1)*y(1))*T234));
  g1(9,8)=1;
  g1(9,14)=(-(T83*getPowerDeriv(y(14),1-params(2),1)));
  g1(9,16)=(-(T86*T234*1/y(1)));
  g1(10,7)=(-1);
  g1(10,9)=1;
  g1(10,10)=(-1);
  g1(10,18)=(-(1+log(y(25))));
  g1(10,25)=(-(y(18)*1/y(25)));
  g1(11,1)=(-(params(1)*y(12)/y(12)*getPowerDeriv(y(1),(-params(4)),1)));
  g1(11,11)=1;
  g1(12,1)=(-((-y(13))*params(8)*(-(y(10)*y(15)))/(y(1)*y(10)*y(1)*y(10))*T245));
  g1(12,10)=(-(T330+(-y(13))*params(8)*T245*(-(y(1)*y(15)))/(y(1)*y(10)*y(1)*y(10))));
  g1(12,12)=1;
  g1(12,13)=(-(T124*(-params(8))));
  g1(12,14)=(-T377);
  g1(12,15)=(-(T330*(-(T111*T113))+(-y(13))*params(8)*T245*1/(y(1)*y(10))));
  g1(13,1)=(-(params(1)*(1-params(8))*(T131*y(13)*getPowerDeriv(y(1),(-params(4)),1)+T129*y(10)/y(15)*T251)));
  g1(13,10)=(-(params(1)*(1-params(8))*T129*T251*y(1)/y(15)+T113*T111*T330));
  g1(13,13)=1-params(1)*(1-params(8))*T104*T131;
  g1(13,14)=(-(T111*T116*getPowerDeriv(y(14),1+params(3),1)+T113*T111*T377));
  g1(13,15)=(-(params(1)*(1-params(8))*T129*T251*(-(y(1)*y(10)))/(y(15)*y(15))+T113*T111*T330*(-(T111*T113))));
  g1(14,9)=(-(T147*1/y(14)));
  g1(14,10)=T142*T330;
  g1(14,12)=T116*y(15)*params(6)*y(14)^params(3)*(-1)/(y(12)*y(12));
  g1(14,14)=T142*T377+T116*1/y(12)*y(15)*params(6)*getPowerDeriv(y(14),params(3),1)-T147*(-y(9))/(y(14)*y(14));
  g1(14,15)=T142*T330*(-(T111*T113))+T116*1/y(12)*params(6)*y(14)^params(3);
  g1(15,1)=(-(y(10)^params(8)*(-y(15))/(y(1)*y(1))*T260));
  g1(15,10)=(-(T153*getPowerDeriv(y(10),params(8),1)));
  g1(15,15)=1-y(10)^params(8)*T260*1/y(1);
  g1(16,1)=(-((1-params(5))*(-y(16))/(y(1)*y(1))));
  g1(16,16)=1-(1-params(5))*1/y(1);
  g1(16,18)=(-1);
  g1(17,1)=(-(log(y(26))-log(y(26))*y(11)*2*T163));
  g1(17,11)=log(y(26))*T163^2;
  g1(17,17)=1;
  g1(17,25)=(-(1/y(25)));
  g1(17,26)=(-(T163*1/y(26)-y(11)*T163^2*1/y(26)));
  g1(18,1)=(-(y(11)*params(2)*y(8)*(params(7)-1)/(params(7)*params(15)*y(16))));
  g1(18,8)=(-(y(11)*params(2)*y(1)*(params(7)-1)/(params(7)*params(15)*y(16))));
  g1(18,11)=(-(params(2)*y(8)*y(1)*(params(7)-1)/(params(7)*params(15)*y(16))+(1-params(5))*y(17)));
  g1(18,16)=(-(y(11)*(-(params(2)*y(8)*y(1)*(params(7)-1)*params(7)*params(15)))/(params(7)*params(15)*y(16)*params(7)*params(15)*y(16))));
  g1(18,17)=1-y(11)*(1-params(5));
  g1(19,19)=1/y(19)-params(12)*1/y(19);
  g1(20,2)=(-y(4));
  g1(20,3)=(-y(4));
  g1(20,4)=(-(y(2)+y(3)-1));
  g1(20,5)=(-1);
  g1(20,16)=(-y(17));
  g1(20,17)=(-y(16));
  g1(20,20)=1;
  g1(20,21)=(-1);
  g1(21,1)=(-(y(11)*(y(21)+y(3)*y(4))));
  g1(21,3)=(-(y(4)*y(1)*y(11)));
  g1(21,4)=(-(y(3)*y(1)*y(11)));
  g1(21,11)=(-(y(1)*(y(21)+y(3)*y(4))));
  g1(21,21)=1-y(1)*y(11);
  g1(22,7)=(-1);
  g1(22,22)=1;
  g1(23,1)=(-(T212*params(16)/2*2*(T67-params(18))));
  g1(23,23)=1;
  g1(24,1)=(-(params(16)*exp(params(16)*(T67-params(18)))));
  g1(24,24)=1;
  g1(25,1)=(-(T222*params(17)/2*2*(T163-params(18))));
  g1(25,25)=1;
  g1(26,1)=(-(params(17)*exp(params(17)*(T163-params(18)))));
  g1(26,26)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],26,676);
end
end
