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

residual = zeros( 24, 1);

%
% Model equations
%

T52 = 1/params(8)*1/params(15);
T59 = (y(2)/y(7))^(1-params(11));
T63 = y(1)*y(7)/y(7);
T78 = (y(16)/y(1))^params(2);
T81 = y(14)^(1-params(2));
T98 = y(1)^(-params(4));
T105 = params(7)/(1+params(3));
T107 = y(14)^(1+params(3));
T109 = y(10)-y(15)*T105*T107;
T110 = T109^(-params(4));
T117 = (y(15)/(y(1)*y(10)))^(1-params(9));
T122 = T98*y(13);
T124 = (y(1)*y(10)/y(15))^params(9);
T135 = y(15)*params(7)*y(14)^params(3)*1/y(12);
T140 = (1-params(2))*1/params(15)*(params(8)-1)/params(8);
T146 = (y(15)/y(1))^(1-params(9));
T156 = y(1)*y(18)/y(18);
T210 = getPowerDeriv(y(16)/y(1),params(2),1);
T221 = getPowerDeriv(y(15)/(y(1)*y(10)),1-params(9),1);
T227 = getPowerDeriv(y(1)*y(10)/y(15),params(9),1);
T236 = getPowerDeriv(y(15)/y(1),1-params(9),1);
T268 = getPowerDeriv(y(2)/y(7),1-params(11),1);
T296 = getPowerDeriv(T109,(-params(4)),1);
T340 = T296*(-(y(15)*T105*getPowerDeriv(y(14),1+params(3),1)));
lhs =y(1);
rhs =params(10)*(1+params(13)*(y(2)+y(3)-1));
residual(1)= lhs-rhs;
lhs =y(1)*y(2);
rhs =params(10)*(y(2)+y(3)*(1-params(13)));
residual(2)= lhs-rhs;
lhs =y(3);
rhs =y(19)*y(2)^(1-params(11))*y(7)^params(11);
residual(3)= lhs-rhs;
lhs =y(4);
rhs =params(13)*y(5)+y(4)*params(10)*(1-params(13))*y(11);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =y(6)+params(10)*y(5)*y(11);
residual(5)= lhs-rhs;
lhs =y(6);
rhs =T52*y(8);
residual(6)= lhs-rhs;
lhs =y(19)*params(11)*y(4)*T59;
rhs =1+y(22)*T63+y(21)-y(11)*y(22)*T63^2;
residual(7)= lhs-rhs;
lhs =y(9);
rhs =y(8);
residual(8)= lhs-rhs;
lhs =y(8);
rhs =T78*T81;
residual(9)= lhs-rhs;
lhs =y(9);
rhs =y(7)+y(10)+(1+y(23))*y(18);
residual(10)= lhs-rhs;
lhs =y(11);
rhs =params(1)*y(12)/y(12)*T98;
residual(11)= lhs-rhs;
lhs =y(12);
rhs =T110+y(13)*params(9)*T117;
residual(12)= lhs-rhs;
lhs =y(13);
rhs =params(1)*(1-params(9))*T122*T124-T107*T105*T110;
residual(13)= lhs-rhs;
lhs =T110*T135;
rhs =T140*y(9)/y(14);
residual(14)= lhs-rhs;
lhs =y(15);
rhs =y(10)^params(9)*T146;
residual(15)= lhs-rhs;
lhs =y(16);
rhs =y(18)+y(16)/y(1)*(1-params(6));
residual(16)= lhs-rhs;
lhs =y(17);
rhs =1+y(23)+T156*y(24)-y(24)*y(11)*T156^2;
residual(17)= lhs-rhs;
lhs =y(17);
rhs =y(11)*(params(2)*y(8)*y(1)*(params(8)-1)/(params(8)*params(15)*y(16))+(1-params(6))*y(17));
residual(18)= lhs-rhs;
lhs =log(y(19));
rhs =log(y(19))*params(14)+0.1*x(1);
residual(19)= lhs-rhs;
lhs =y(20);
rhs =y(3)*y(4);
residual(20)= lhs-rhs;
lhs =y(21);
rhs =params(16)/2*(T63-params(18))^2;
residual(21)= lhs-rhs;
lhs =y(22);
rhs =params(16)*(T63-params(18));
residual(22)= lhs-rhs;
lhs =y(23);
rhs =params(17)/2*(T156-params(18))^2;
residual(23)= lhs-rhs;
lhs =y(24);
rhs =params(17)*(T156-params(18));
residual(24)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(24, 24);

  %
  % Jacobian matrix
  %

  g1(1,1)=1;
  g1(1,2)=(-(params(10)*params(13)));
  g1(1,3)=(-(params(10)*params(13)));
  g1(2,1)=y(2);
  g1(2,2)=y(1)-params(10);
  g1(2,3)=(-(params(10)*(1-params(13))));
  g1(3,2)=(-(y(7)^params(11)*y(19)*getPowerDeriv(y(2),1-params(11),1)));
  g1(3,3)=1;
  g1(3,7)=(-(y(19)*y(2)^(1-params(11))*getPowerDeriv(y(7),params(11),1)));
  g1(3,19)=(-(y(2)^(1-params(11))*y(7)^params(11)));
  g1(4,4)=1-params(10)*(1-params(13))*y(11);
  g1(4,5)=(-params(13));
  g1(4,11)=(-(y(4)*params(10)*(1-params(13))));
  g1(5,5)=1-params(10)*y(11);
  g1(5,6)=(-1);
  g1(5,11)=(-(params(10)*y(5)));
  g1(6,6)=1;
  g1(6,8)=(-T52);
  g1(7,1)=(-(y(22)-y(11)*y(22)*2*T63));
  g1(7,2)=y(19)*params(11)*y(4)*1/y(7)*T268;
  g1(7,4)=T59*y(19)*params(11);
  g1(7,7)=y(19)*params(11)*y(4)*T268*(-y(2))/(y(7)*y(7));
  g1(7,11)=y(22)*T63^2;
  g1(7,19)=params(11)*y(4)*T59;
  g1(7,21)=(-1);
  g1(7,22)=(-(T63-y(11)*T63^2));
  g1(8,8)=(-1);
  g1(8,9)=1;
  g1(9,1)=(-(T81*(-y(16))/(y(1)*y(1))*T210));
  g1(9,8)=1;
  g1(9,14)=(-(T78*getPowerDeriv(y(14),1-params(2),1)));
  g1(9,16)=(-(T81*T210*1/y(1)));
  g1(10,7)=(-1);
  g1(10,9)=1;
  g1(10,10)=(-1);
  g1(10,18)=(-(1+y(23)));
  g1(10,23)=(-y(18));
  g1(11,1)=(-(params(1)*y(12)/y(12)*getPowerDeriv(y(1),(-params(4)),1)));
  g1(11,11)=1;
  g1(12,1)=(-(y(13)*params(9)*(-(y(10)*y(15)))/(y(1)*y(10)*y(1)*y(10))*T221));
  g1(12,10)=(-(T296+y(13)*params(9)*T221*(-(y(1)*y(15)))/(y(1)*y(10)*y(1)*y(10))));
  g1(12,12)=1;
  g1(12,13)=(-(params(9)*T117));
  g1(12,14)=(-T340);
  g1(12,15)=(-(T296*(-(T105*T107))+y(13)*params(9)*T221*1/(y(1)*y(10))));
  g1(13,1)=(-(params(1)*(1-params(9))*(T124*y(13)*getPowerDeriv(y(1),(-params(4)),1)+T122*y(10)/y(15)*T227)));
  g1(13,10)=(-(params(1)*(1-params(9))*T122*T227*y(1)/y(15)-T107*T105*T296));
  g1(13,13)=1-params(1)*(1-params(9))*T98*T124;
  g1(13,14)=T105*T110*getPowerDeriv(y(14),1+params(3),1)+T107*T105*T340;
  g1(13,15)=(-(params(1)*(1-params(9))*T122*T227*(-(y(1)*y(10)))/(y(15)*y(15))-T107*T105*T296*(-(T105*T107))));
  g1(14,9)=(-(T140*1/y(14)));
  g1(14,10)=T135*T296;
  g1(14,12)=T110*y(15)*params(7)*y(14)^params(3)*(-1)/(y(12)*y(12));
  g1(14,14)=T135*T340+T110*1/y(12)*y(15)*params(7)*getPowerDeriv(y(14),params(3),1)-T140*(-y(9))/(y(14)*y(14));
  g1(14,15)=T135*T296*(-(T105*T107))+T110*1/y(12)*params(7)*y(14)^params(3);
  g1(15,1)=(-(y(10)^params(9)*(-y(15))/(y(1)*y(1))*T236));
  g1(15,10)=(-(T146*getPowerDeriv(y(10),params(9),1)));
  g1(15,15)=1-y(10)^params(9)*T236*1/y(1);
  g1(16,1)=(-((1-params(6))*(-y(16))/(y(1)*y(1))));
  g1(16,16)=1-(1-params(6))*1/y(1);
  g1(16,18)=(-1);
  g1(17,1)=(-(y(24)-y(24)*y(11)*2*T156));
  g1(17,11)=y(24)*T156^2;
  g1(17,17)=1;
  g1(17,23)=(-1);
  g1(17,24)=(-(T156-y(11)*T156^2));
  g1(18,1)=(-(y(11)*params(2)*y(8)*(params(8)-1)/(params(8)*params(15)*y(16))));
  g1(18,8)=(-(y(11)*params(2)*y(1)*(params(8)-1)/(params(8)*params(15)*y(16))));
  g1(18,11)=(-(params(2)*y(8)*y(1)*(params(8)-1)/(params(8)*params(15)*y(16))+(1-params(6))*y(17)));
  g1(18,16)=(-(y(11)*(-(params(2)*y(8)*y(1)*(params(8)-1)*params(8)*params(15)))/(params(8)*params(15)*y(16)*params(8)*params(15)*y(16))));
  g1(18,17)=1-y(11)*(1-params(6));
  g1(19,19)=1/y(19)-params(14)*1/y(19);
  g1(20,3)=(-y(4));
  g1(20,4)=(-y(3));
  g1(20,20)=1;
  g1(21,1)=(-(params(16)/2*2*(T63-params(18))));
  g1(21,21)=1;
  g1(22,1)=(-params(16));
  g1(22,22)=1;
  g1(23,1)=(-(params(17)/2*2*(T156-params(18))));
  g1(23,23)=1;
  g1(24,1)=(-params(17));
  g1(24,24)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],24,576);
end
end
