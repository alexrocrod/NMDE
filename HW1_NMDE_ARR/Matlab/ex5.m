%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021

close all
clear all
%% Question 5
% Lotka–Volterra

% x'(t) = x(t) (α − βy(t))
% y'(t) = y(t) (γx(t) − δ)
% x(0) = x0
% y(0) = y0

%% Constants

alpha = 0.2;
beta = 0.01;
gamma = 0.004;
delta = 0.07;

x0 = 19;
y0 = 22;

%% Variables

t0 = 0;
T = 300;
h = 10^-3;
ts = t0:h:T;
N = length(ts);

f = @(r) [r(1) * (alpha - beta*r(2))     r(2) * (gamma*r(1) - delta)];

r = zeros(N,2);
r(1,1) = x0;
r(1,2) = y0;

%% Solver

tic
for i=2:N
    k1 = f(r(i-1,:));
    k2 = f(r(i-1,:)+h/2*k1);
    k3 = f(r(i-1,:)+h/2*k2);
    k4 = f(r(i-1,:)+h*k3);

    r(i,:) = r(i-1,:) + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
end
toc

%% Plot

plot(ts, r(:,1), 'b.', ts, r(:,2), 'r.')
legend("prey","predator",'Location','best')

