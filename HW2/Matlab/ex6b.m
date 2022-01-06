%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 6c) - Solving 3b)

n = 1e4;
v = ones(n,1);
vi = 1:5;
v(vi) = 200*vi;

A = sparse(diag(v));
L = ichol(A);
n = size(A, 1);
b = rand(n, 1);
tol = 1e-8;
maxit = 200;
restart = maxit*10;

tic
% [x, flag, relres, iter1, resvec1] = gmres( A, b, restart, tol, maxit,L,L');
[x, flag, relres, iter1, resvec1] = gmres( A, b, restart, tol, maxit)
toc
totalit = (iter1(1)-1)*restart + iter1(2);

L = speye(size(L));
tic
[x, resvec, iter] = mypcg(A, b, tol, maxit, L);
toc

semilogy(0:totalit, resvec1, 'r-*',0:iter, resvec, 'g-+')
legend('Matlab GMRES' , 'My PCG');
xlabel('Iterations');
ylabel('Residual Norm');

