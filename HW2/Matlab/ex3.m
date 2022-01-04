%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 3
n = 1e4;
v = ones(n,1);
vi = 1:5;
v(vi) = 200*vi;

A = sparse(diag(v));
k = condest(A);
L = ichol(A);
n = size(A, 1);
b = rand(n, 1);
tol = 1e-10;
maxit = 200;

% tic
% [x, flag, relres, iter1, resvec1] = pcg( A, b, tol, maxit,L,L')
% toc
tic
[x, flag, relres, iter1, resvec1] = pcg( A, b, tol, maxit);
toc

L = eye(size(L));
tic
[x, resvec, iter] = mypcg(A, b, tol, maxit, L);
toc

semilogy(0:iter1, resvec1, 'g-o',0:iter, resvec, 'r-*')
xlabel('Iterations');
ylabel('Residual Norm');

