%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021

close all
clear all

%% Question 1
% 2-step Simpson's method:
% y(n+2) = y(n) + h/3 * (f(n) + 4*f(n+1) + f(n+2));
% y' = -5y;

%% Constants
h = 0.02;
T = 10;
ts = 0:h:T;

N = length(ts);
y = zeros(N,1);
y(1) = 1;
f = @(y) -5*y; 
y_exact = exp(-5.*ts);

%% Compute y(2) - FE
y(2) = y(1) + h*f(y(1));

%% Compute y(2) - RK4
% k1 = f(y(1));
% k2 = f(y(1)+h/2*k1);
% k3 = f(y(1)+h/2*k2);
% k4 = f(y(1)+h*k3);
% 
% y(2) = y(1) + h/6 * (k1 + 2*k2 + 2*k3 + k4);

%% Main Method
for i=2:N-2
    y(i+2) = y(i) + h/3 * (f(y(i)) + 4*f(y(i+1)) + f(y(i+2)));
end

% RK4 for extra comparison
% for i=3:N
%         k1 = f(y(i-1));
%         k2 = f(y(i-1)+h/2*k1);
%         k3 = f(y(i-1)+h/2*k2);
%         k4 = f(y(i-1)+h*k3);
%         
%         y(i) = y(i-1) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
% end

y=y';
error = abs(y-y_exact);
plot(ts,error,'.')
xlabel("time")
ylabel("absolute error")

errorf = error(end)

