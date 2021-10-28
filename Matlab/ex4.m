%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021

close all
clear all

%% Question 4

% y' = -Ay;
%% Constants
nx = 100;
G = numgrid ( 'S' , nx ) ;
A = delsq(G) * ( nx-1)^2 ;
lambda = -eigs(A,1,'lm');

h=lambda^-3;
T = 0.1;
ts = h:h:T;

N = length(ts);
y = zeros(N,N);
y(1) = ones(1,N);

%y_exact_f = exp(0.1*A)*y(1);

%% Main

for i=2:N
    k1 = f(y(i-1));
    k2 = f(y(i-1)+h/2*k1);
    k3 = f(y(i-1)+h/2*k2);
    k4 = f(y(i-1)+h*k3);
        
    y(i) = y(i-1) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

%% 



