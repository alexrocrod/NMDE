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
% y'(t) = y(t) (δ − γx(t))
% x(0) = x0
% y(0) = y0

%% Constants
alpha = 0.2;
beta = 0.01;
gama = 0.07;
delta = 0.004;

x0 = 19;
y0 = 22;
t0 = 0;
T = 300;
h = 10^-3;
ts = t0:h:T;
N = length(ts);

fx = @(x,y) x * (alpha - beta*y);
fy = @(x,y) y * (delta - gama*x);

x = zeros(N,1);
x(1) = x0;
y = zeros(N,1);
y(1) = y0;


%% Main 

for i=2:N
    %x
%     k1 = fx(x(i-1),y(i-1));
%     k2 = fx(x(i-1)+h/2*k1,y(i-1)+h/2*k1);
%     k3 = fx(x(i-1)+h/2*k2,y(i-1)+h/2*k2);
%     k4 = fx(x(i-1)+h*k3,y(i-1)+h*k3);

    k1 = fx(x(i-1),y(i-1));
    k2 = fx(x(i-1)+h/2*k1,y(i-1));
    k3 = fx(x(i-1)+h/2*k2,y(i-1));
    k4 = fx(x(i-1)+h*k3,y(i-1));
        
    x(i) = x(i-1) + h/6 * (k1 + 2*k2 + 2*k3 + k4);

    %y

%     k1 = fy(x(i-1),y(i-1));
%     k2 = fy(x(i-1)+h/2*k1,y(i-1)+h/2*k1);
%     k3 = fy(x(i-1)+h/2*k2,y(i-1)+h/2*k2);
%     k4 = fy(x(i-1)+h*k3,y(i-1)+h*k3);

    k1 = fy(x(i-1),y(i-1));
    k2 = fy(x(i-1),y(i-1)+h/2*k1);
    k3 = fy(x(i-1),y(i-1)+h/2*k2);
    k4 = fy(x(i-1),y(i-1)+h*k3);
        
    y(i) = y(i-1) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

plot(ts,x,'b.')
hold on
plot(ts,y,'r.')









