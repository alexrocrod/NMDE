%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 1
% Implementation of the PCG method with Cholesky preconditioner

A = delsq(numgrid('S', 102 ));
L = ichol(A);
n = size(A, 1);
b = A * ones(n, 1);
tol = 1e-8;
maxit = 250; % 750 to understand convergence of my PCG implementation

%% Main

% Matlab PCG Without Preconditioner
tic
[x1, flag1, relres1, iter1, resvec1] = pcg(A, b, tol, maxit); 
toc

% Matlab PCG With IC(0)
tic
[x2, flag2, relres2, iter2, resvec2] = pcg(A, b, tol, maxit, L, L');
toc

% My PCG Implementation
tic
[x3, resvec3, iter3] = mypcg(A, b, tol, maxit, L);
toc

semilogy(0:iter1, resvec1, 'r-*', 0:iter2, resvec2,'g-o', 0:iter3, resvec3,'b-+')
legend('No preconditioner' , 'IC(0)', 'My implementation');
xlabel('Iterations');
ylabel('Residual Norm');


