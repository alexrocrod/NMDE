%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021

close all
clear all


load("exact_solution.txt")
t_exact = exact_solution(:,1);
y_exact = exact_solution(:,2);


%% Question 4

% y' = -Ay;
%% Constants
nx = 100;
G = numgrid ( 'S' , nx ) ;
A = delsq(G) * ( nx-1)^2 ;
lambda = -eigs(A,1,'lm');

% h=(-lambda)^-1;
T = 0.1;
h=T/length(A);
ts = h:h:T;

N = length(ts);
y = zeros(N,N);
y(1,:) = ones(1,N);

%y_exact_f = exp(0.1*A)*y(1);
f = @(x) -A*x';

%% Stability with RK4
% tic
% for i=2:N
%     k1 = f(y(i-1,:));
%     k2 = f(y(i-1,:)+k1'.*h/2);
%     k3 = f(y(i-1,:)+h/2.*k2');
%     k4 = f(y(i-1,:)+h.*k3');
%         
%     y(i,:) = y(i-1,:) + h/6 * (k1 + 2*k2 + 2*k3 + k4)';
% end
% toc

%% Solve with ODE45
tic
tspan = [0 0.1];
[t,y] = ode45(@(t,y) -A*y, tspan, y(1,:));
toc

yf = y(end,:)';
% plot(ts, yf,'.')
% figure(2)
% plot(ts, y_exact,'.')
figure(3)
plot(ts, abs(yf-y_exact))


%% Solve with CN
ords=3:5;
hs1=10.^-ords;
for h=hs1
    tol = h^-3;



end


%% Solve with BDF3





