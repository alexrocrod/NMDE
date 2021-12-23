%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 1
% Implementation of the PCG method with Cholesky preconditioner

%% Test Variables
A = delsq(numgrid('S', 102 ));
L = ichol(A);
n = size(A, 1);
b=A * ones(n, 1);
tol = 1e-8;
maxit = 200;

%% Main

% Matlab PCG without preconditioner
tic
[x, flag, relres, iter1, resvec1] = pcg( A, b, tol, maxit); 
toc
% Matlab PCG with IC(0)
tic
[x, flag, relres, iter2, resvec2] = pcg( A, b, tol, maxit, L, L');
toc

% my implementation
tic
[x, resvec3, iter3] = mypcg(A, b, tol, maxit, L);
toc
semilogy(0:iter1, resvec1, 'r-*', 0:iter2, resvec2,'g-o', 0:iter3, resvec3,'b-+')
legend('No preconditioner' , 'IC(0)', 'My implementation');
xlabel('Iterations');
ylabel('Residual Norm');


