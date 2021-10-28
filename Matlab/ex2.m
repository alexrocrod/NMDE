%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021

close all
clear all

%% Question 2
% 4-stage Runge-Kutta method:
% y(n+2) = y(n) + h/3 * (f(n) + 4*f(n+1) + f(n+2));

% y' = -10y^2;
%% Constants
ks= 5:10;
T = 2;
hs = 2.^-ks;
errorf=zeros(length(hs),1);
Nsteps=zeros(length(hs),1);
ie=1;
for h=hs
    ts = 0:h:T;
    
    N = length(ts);
    y = zeros(N,1);
    y(1) = 1;
    
    f = @(y) -10*y^2; 
    y_exact = -10/3*ts.^3 + 1;
    
    % Main Method
    for i=2:N
        k1 = f(y(i-1));
        k2 = f(y(i-1)+h/2*k1);
        k3 = f(y(i-1)+h/2*k2);
        k4 = f(y(i-1)+h*k3);
        
        y(i) = y(i-1) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
    
    Nsteps(ie) = N;
    errorf(ie) = abs(y(end)-y_exact(T));
    ie=ie+1;
   
end
   
plot(log(Nsteps),log(errorf))

tab = [transpose(hs) errorf]

