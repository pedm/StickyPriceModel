%% Define parameter bounds
% Search over 2^9 grid of parameters

LB.eta = 0;  
LB.gamma = 0.00001;
LB.phi = 0.5;
LB.lambda_bar = 0.005;
LB.psi_N = 0;
LB.rhozeta = 0.0001; 
LB.rhozeta2 = 0.0001;  
LB.sigmazeta = .5;
LB.zetabar = .1;
LB.rho_lambda = .01;

UB.eta = 0.9999;  
UB.gamma = 0.9999; 
UB.phi = 0.99; 
UB.lambda_bar = 10; 
UB.psi_N = 100; 
UB.rhozeta = 1.5; 
UB.rhozeta2 = 0.99; 
UB.sigmazeta = 10;  
UB.zetabar = 10;
UB.rho_lambda = .99;

%% Not actual bounds, but useful restrictions
LB.rhozeta2 = 0.0001;
UB.rhozeta2 = 0.2;

%% Set distance between

FF = fields(LB);
 
DIST = [];
DIST.(char(FF(1))) = UB.(char(FF(1))) - LB.(char(FF(1)));
DIST.(char(FF(2))) = UB.(char(FF(2))) - LB.(char(FF(2)));
DIST.(char(FF(3))) = UB.(char(FF(3))) - LB.(char(FF(3)));
DIST.(char(FF(4))) = UB.(char(FF(4))) - LB.(char(FF(4)));
DIST.(char(FF(5))) = UB.(char(FF(5))) - LB.(char(FF(5)));
DIST.(char(FF(6))) = UB.(char(FF(6))) - LB.(char(FF(6)));
DIST.(char(FF(7))) = UB.(char(FF(7))) - LB.(char(FF(7)));
DIST.(char(FF(8))) = UB.(char(FF(8))) - LB.(char(FF(8)));
DIST.(char(FF(9))) = UB.(char(FF(9))) - LB.(char(FF(9)));
DIST.(char(FF(10))) = UB.(char(FF(10))) - LB.(char(FF(10)));


guess = [];
guess(1) =  LB.(char(FF(1))) + DIST.(char(FF(1)))*rand(1,1);
guess(2) =  LB.(char(FF(2))) + DIST.(char(FF(2)))*rand(1,1);
guess(3) =  LB.(char(FF(3))) + DIST.(char(FF(3)))*rand(1,1);
guess(4) =  LB.(char(FF(4))) + DIST.(char(FF(4)))*rand(1,1);
guess(5) =  LB.(char(FF(5))) + DIST.(char(FF(5)))*rand(1,1);
guess(6) =  LB.(char(FF(6))) + DIST.(char(FF(6)))*rand(1,1);
guess(7) =  LB.(char(FF(7))) + DIST.(char(FF(7)))*rand(1,1);
guess(8) =  LB.(char(FF(8))) + DIST.(char(FF(8)))*rand(1,1);
guess(9) =  LB.(char(FF(10))) + DIST.(char(FF(10)))*rand(1,1); % zetabar is not being guessed, but rho_lambda is
