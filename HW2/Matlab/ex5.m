%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 5
A = load("mat13041.rig");
A = spconvert(A);
n = size(A, 1);
b1 = 1./sqrt(1:n);
x_exact = b1';
b = A * x_exact;
tol = 1e-10;
maxit = 550;
x0 = ones(n,1);

% Matlab GMRES without restart
restart = n;
tic
[~, ~, ~, iter1, resvec1] = gmres(A, b, restart, tol, maxit);
toc
totalit = (iter1(1)-1)*restart + iter1(2);

% my GMRES
tic
[~, iter2, resvec2, ~] = mygmres(A, b, tol, maxit, x0);
toc

semilogy(0:totalit, resvec1, 'r-*', 0:iter2, resvec2', 'g-+')
legend('Matlab GMRES' , 'My GMRES');
xlabel('Iterations');
ylabel('Residual Norm');

