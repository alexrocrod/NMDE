%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021


%% Question 4

% y' = -Ay;
%% Constants
nx = 100;
G = numgrid ( 'S' , nx ) ;
A = delsq(G) * ( nx-1)^2 ;

h=?;
T = 0.1;
ts = 0:h:T;

N = length(ts);
y = Ones(N+1,1);
y(1) = 1;

lambda = -eigs(A,1,'lm');

%% 