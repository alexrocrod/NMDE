%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021


%% Question 1
% 2-step Simpson's method:
% y(n+2) = y(n) + h/3 * (f(n) + 4*f(n+1) + f(n+2));

% y' = -5y;
%% Constants
h = 0.02;
T = 10;
ts = 0:h:T;

N = length(ts);
y = Zeros(N+1,1);
y(1) = 1;



%% 