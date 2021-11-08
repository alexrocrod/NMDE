%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021

close all
clear all

%% Question 1 e)
% find the optimal y(1) = alfa to minimize the error

% 2-step Simpson's method:
% y(n+2) = y(n) + h/3 * (f(n) + 4*f(n+1) + f(n+2));
% y' = -5y;

%% Constants
h = 0.02;
T = 10;
ts = 0:h:T;

N = length(ts);
y = zeros(N,1);

f = @(y) -5*y; % function for ys not in the vector 
y_exact = -5/2*ts.^2; % without constant (alfa)

ial=1;
% alfas=0.0001:0.0001:1; % finds that it is around 0.0025
% alfas=0.002:0.000001:0.003; % finds a more precise value
alfas=0.00200:0.00000001:0.00205; % finds a even more precise value

y_al = zeros(1,length(alfas));
for alfa=alfas
    y(1) = alfa;
    %% Compute y(2) - RK4
    k1 = f(y(1));
    k2 = f(y(1)+h/2*k1);
    k3 = f(y(1)+h/2*k2);
    k4 = f(y(1)+h*k3);
    
    y(2) = y(1) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    %% Main Method
    for i=1:N-2
        y(i+2) = y(i) + h/3 * (f(y(i)) + 4*f(y(i+1)) + f(y(i+2)));
    end
    y_al(ial) = y(end);
    ial = ial+1;
end

y_f_ex = ones(1,length(alfas)) * y_exact(end) + alfas;

plot(alfas,abs(y_al-y_f_ex),'.')

[M,I] = min(abs(y_al-y_f_ex))

alfaf = alfas(I);

format long
fprintf("minimum error = %d , alfa = %d\n", M, alfaf)

figure(2)

plot(alfas,abs(y_al-y_f_ex),'.')
ylim([0 10*M])
xlim([alfas(1) alfas(end)])



