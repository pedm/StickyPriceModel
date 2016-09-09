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
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 26, 1);

%
% Model equations
%

T52 = 1/params(8)*1/params(14);
T59 = (y(2)/y(7))^(1-params(11));
T64 = y(1)*y(7)/y(7);
T80 = (y(16)/y(1))^params(2);
T83 = y(14)^(1-params(2));
T101 = y(1)^(-params(4));
T108 = params(7)/(1+params(3));
T110 = y(14)^(1+params(3));
T112 = y(10)-y(15)*T108*T110;
T113 = T112^(-params(4));
T121 = (y(15)/(y(1)*y(10)))^(1-params(9));
T126 = T101*y(13);
T128 = (y(1)*y(10)/y(15))^params(9);
T139 = y(15)*params(7)*y(14)^params(3)*1/y(12);
T144 = (1-params(2))*1/params(14)*(params(8)-1)/params(8);
T150 = (y(15)/y(1))^(1-params(9));
T160 = y(1)*y(18)/y(18);
T209 = exp(params(15)/2*(T64-params(17))^2);
T219 = exp(params(16)/2*(T160-params(17))^2);
T231 = getPowerDeriv(y(16)/y(1),params(2),1);
T242 = getPowerDeriv(y(15)/(y(1)*y(10)),1-params(9),1);
T248 = getPowerDeriv(y(1)*y(10)/y(15),params(9),1);
T257 = getPowerDeriv(y(15)/y(1),1-params(9),1);
T295 = getPowerDeriv(y(2)/y(7),1-params(11),1);
T328 = getPowerDeriv(T112,(-params(4)),1);
T375 = T328*(-(y(15)*T108*getPowerDeriv(y(14),1+params(3),1)));
lhs =y(1);
rhs =params(10)*(1+params(12)*(y(2)+y(3)-1));
residual(1)= lhs-rhs;
lhs =y(1)*y(2);
rhs =params(10)*(y(2)+y(3)*(1-params(12)));
residual(2)= lhs-rhs;
lhs =y(3);
rhs =y(19)*y(2)^(1-params(11))*y(7)^params(11);
residual(3)= lhs-rhs;
lhs =y(4);
rhs =params(12)*y(5)+y(4)*params(10)*(1-params(12))*y(11);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =y(6)+params(10)*y(5)*y(11);
residual(5)= lhs-rhs;
lhs =y(6);
rhs =T52*y(8);
residual(6)= lhs-rhs;
lhs =y(19)*params(11)*y(4)*T59;
rhs =1+log(y(24))*T64+log(y(23))-y(11)*log(y(24))*T64^2;
residual(7)= lhs-rhs;
lhs =y(9);
rhs =y(8);
residual(8)= lhs-rhs;
lhs =y(8);
rhs =T80*T83;
residual(9)= lhs-rhs;
lhs =y(9);
rhs =y(7)+y(10)+(1+log(y(25)))*y(18);
residual(10)= lhs-rhs;
lhs =y(11);
rhs =params(1)*y(12)/y(12)*T101;
residual(11)= lhs-rhs;
lhs =y(12);
rhs =T113+(-y(13))*params(9)*T121;
residual(12)= lhs-rhs;
lhs =y(13);
rhs =params(1)*(1-params(9))*T126*T128+T110*T108*T113;
residual(13)= lhs-rhs;
lhs =T113*T139;
rhs =T144*y(9)/y(14);
residual(14)= lhs-rhs;
lhs =y(15);
rhs =y(10)^params(9)*T150;
residual(15)= lhs-rhs;
lhs =y(16);
rhs =y(18)+y(16)/y(1)*(1-params(6));
residual(16)= lhs-rhs;
lhs =y(17);
rhs =1+log(y(25))+T160*log(y(26))-log(y(26))*y(11)*T160^2;
residual(17)= lhs-rhs;
lhs =y(17);
rhs =y(11)*(params(2)*y(8)*y(1)*(params(8)-1)/(params(8)*params(14)*y(16))+(1-params(6))*y(17));
residual(18)= lhs-rhs;
lhs =log(y(19));
rhs =log(y(19))*params(13)+0.1*x(1);
residual(19)= lhs-rhs;
lhs =y(20);
rhs =y(5)+y(16)*y(17)+(y(2)+y(3)-1)*y(4)+y(21);
residual(20)= lhs-rhs;
lhs =y(21);
rhs =y(1)*y(11)*(y(21)+y(3)*y(4));
residual(21)= lhs-rhs;
lhs =y(22);
rhs =y(3)*y(4);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =T209;
residual(23)= lhs-rhs;
lhs =y(24);
rhs =exp(params(15)*(T64-params(17)));
residual(24)= lhs-rhs;
lhs =y(25);
rhs =T219;
residual(25)= lhs-rhs;
lhs =y(26);
rhs =exp(params(16)*(T160-params(17)));
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
  g1(1,2)=(-(params(10)*params(12)));
  g1(1,3)=(-(params(10)*params(12)));
  g1(2,1)=y(2);
  g1(2,2)=y(1)-params(10);
  g1(2,3)=(-(params(10)*(1-params(12))));
  g1(3,2)=(-(y(7)^params(11)*y(19)*getPowerDeriv(y(2),1-params(11),1)));
  g1(3,3)=1;
  g1(3,7)=(-(y(19)*y(2)^(1-params(11))*getPowerDeriv(y(7),params(11),1)));
  g1(3,19)=(-(y(2)^(1-params(11))*y(7)^params(11)));
  g1(4,4)=1-params(10)*(1-params(12))*y(11);
  g1(4,5)=(-params(12));
  g1(4,11)=(-(y(4)*params(10)*(1-params(12))));
  g1(5,5)=1-params(10)*y(11);
  g1(5,6)=(-1);
  g1(5,11)=(-(params(10)*y(5)));
  g1(6,6)=1;
  g1(6,8)=(-T52);
  g1(7,1)=(-(log(y(24))-y(11)*log(y(24))*2*T64));
  g1(7,2)=y(19)*params(11)*y(4)*1/y(7)*T295;
  g1(7,4)=T59*y(19)*params(11);
  g1(7,7)=y(19)*params(11)*y(4)*T295*(-y(2))/(y(7)*y(7));
  g1(7,11)=log(y(24))*T64^2;
  g1(7,19)=params(11)*y(4)*T59;
  g1(7,23)=(-(1/y(23)));
  g1(7,24)=(-(T64*1/y(24)-T64^2*y(11)*1/y(24)));
  g1(8,8)=(-1);
  g1(8,9)=1;
  g1(9,1)=(-(T83*(-y(16))/(y(1)*y(1))*T231));
  g1(9,8)=1;
  g1(9,14)=(-(T80*getPowerDeriv(y(14),1-params(2),1)));
  g1(9,16)=(-(T83*T231*1/y(1)));
  g1(10,7)=(-1);
  g1(10,9)=1;
  g1(10,10)=(-1);
  g1(10,18)=(-(1+log(y(25))));
  g1(10,25)=(-(y(18)*1/y(25)));
  g1(11,1)=(-(params(1)*y(12)/y(12)*getPowerDeriv(y(1),(-params(4)),1)));
  g1(11,11)=1;
  g1(12,1)=(-((-y(13))*params(9)*(-(y(10)*y(15)))/(y(1)*y(10)*y(1)*y(10))*T242));
  g1(12,10)=(-(T328+(-y(13))*params(9)*T242*(-(y(1)*y(15)))/(y(1)*y(10)*y(1)*y(10))));
  g1(12,12)=1;
  g1(12,13)=(-(T121*(-params(9))));
  g1(12,14)=(-T375);
  g1(12,15)=(-(T328*(-(T108*T110))+(-y(13))*params(9)*T242*1/(y(1)*y(10))));
  g1(13,1)=(-(params(1)*(1-params(9))*(T128*y(13)*getPowerDeriv(y(1),(-params(4)),1)+T126*y(10)/y(15)*T248)));
  g1(13,10)=(-(params(1)*(1-params(9))*T126*T248*y(1)/y(15)+T110*T108*T328));
  g1(13,13)=1-params(1)*(1-params(9))*T101*T128;
  g1(13,14)=(-(T108*T113*getPowerDeriv(y(14),1+params(3),1)+T110*T108*T375));
  g1(13,15)=(-(params(1)*(1-params(9))*T126*T248*(-(y(1)*y(10)))/(y(15)*y(15))+T110*T108*T328*(-(T108*T110))));
  g1(14,9)=(-(T144*1/y(14)));
  g1(14,10)=T139*T328;
  g1(14,12)=T113*y(15)*params(7)*y(14)^params(3)*(-1)/(y(12)*y(12));
  g1(14,14)=T139*T375+T113*1/y(12)*y(15)*params(7)*getPowerDeriv(y(14),params(3),1)-T144*(-y(9))/(y(14)*y(14));
  g1(14,15)=T139*T328*(-(T108*T110))+T113*1/y(12)*params(7)*y(14)^params(3);
  g1(15,1)=(-(y(10)^params(9)*(-y(15))/(y(1)*y(1))*T257));
  g1(15,10)=(-(T150*getPowerDeriv(y(10),params(9),1)));
  g1(15,15)=1-y(10)^params(9)*T257*1/y(1);
  g1(16,1)=(-((1-params(6))*(-y(16))/(y(1)*y(1))));
  g1(16,16)=1-(1-params(6))*1/y(1);
  g1(16,18)=(-1);
  g1(17,1)=(-(log(y(26))-log(y(26))*y(11)*2*T160));
  g1(17,11)=log(y(26))*T160^2;
  g1(17,17)=1;
  g1(17,25)=(-(1/y(25)));
  g1(17,26)=(-(T160*1/y(26)-y(11)*T160^2*1/y(26)));
  g1(18,1)=(-(y(11)*params(2)*y(8)*(params(8)-1)/(params(8)*params(14)*y(16))));
  g1(18,8)=(-(y(11)*params(2)*y(1)*(params(8)-1)/(params(8)*params(14)*y(16))));
  g1(18,11)=(-(params(2)*y(8)*y(1)*(params(8)-1)/(params(8)*params(14)*y(16))+(1-params(6))*y(17)));
  g1(18,16)=(-(y(11)*(-(params(2)*y(8)*y(1)*(params(8)-1)*params(8)*params(14)))/(params(8)*params(14)*y(16)*params(8)*params(14)*y(16))));
  g1(18,17)=1-y(11)*(1-params(6));
  g1(19,19)=1/y(19)-params(13)*1/y(19);
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
  g1(22,3)=(-y(4));
  g1(22,4)=(-y(3));
  g1(22,22)=1;
  g1(23,1)=(-(T209*params(15)/2*2*(T64-params(17))));
  g1(23,23)=1;
  g1(24,1)=(-(params(15)*exp(params(15)*(T64-params(17)))));
  g1(24,24)=1;
  g1(25,1)=(-(T219*params(16)/2*2*(T160-params(17))));
  g1(25,25)=1;
  g1(26,1)=(-(params(16)*exp(params(16)*(T160-params(17)))));
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
