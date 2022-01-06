%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 6
% 
load("mat13041.rig");
A = spconvert(mat13041);
n = size(A, 1);
b1 = 1./sqrt(1:n);
x_exact = b1';
b = A * x_exact;
tol = 1e-10; %1e-10;
maxit = 550;
x0 = zeros(n,1);
restart = n;

setup.type = 'crout';
setup.droptol = 0.01;
tic
% [L,U] = ilu(A);
[L,U] = ilu(A,setup);
toc

tic
[x1, flag1, relres, iter1, resvec1] = gmres(A, b, restart, tol, maxit,L,U);
toc
totalit = (iter1(1)-1)*restart + iter1(2);
norm(b-A*x1)

tic
[x2, iter2, resvec2, flag2] = myprecgmres(A, b, tol, maxit, x0, L, U);
toc
norm(b-A*x2)

semilogy(0:totalit, resvec1, 'r-*',0:iter2, resvec2, 'g-+')
legend('Matlab GMRES' , 'My Prec. GMRES');
xlabel('Iterations');
ylabel('Residual Norm');



