%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 4
A = gallery('wathen',100,100);
n = size(A, 1);
b = rand(n,1);
tol = 1e-8;
maxit = 500;


% non preconditioned
tic
[x, flag, relres, iter1, resvec1] = pcg( A, b, tol, maxit); 
toc

% Jacobi
M = sparse(diag(diag(A)));
tic
[x, flag, relres, iter2, resvec2] = pcg( A, b, tol, maxit, M); 
toc

% IC(0)
L = ichol(A);
tic
[x, flag, relres, iter3, resvec3] = pcg( A, b, tol, maxit, L, L'); 
toc

semilogy(0:iter1, resvec1, 'r-*', 0:iter2, resvec2,'g-o', 0:iter3, resvec3,'b-+')
legend('No preconditioner', 'Jacobi', 'IC(0)');
xlabel('Iterations');
ylabel('Residual Norm');

