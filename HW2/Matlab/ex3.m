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
v(vi) = 200 * vi;

A = sparse(diag(v));
k = condest(A);
L = ichol(A);
n = size(A, 1);
x_exact = rand(n,1);
b = A * x_exact;
tol = 1e-8;
maxit = 20;

%% Original Exercise (only 1 iteration in both methods)

tic
[~,~,~, iter1, resvec1] = pcg(A, b, tol, maxit,L,L');
toc

tic
[~, resvec2, iter2] = mypcg(A, b, tol, maxit, L);
toc

semilogy(0:iter1, resvec1, 'g-o', 0:iter2, resvec2, 'r-*')
legend('Matlab PCG IC(0)', 'My PCG IC(0)');
xlabel('Iterations');
ylabel('Residual Norm');


%% Without Preconditioning (larger number of iterations)

tic
[~,~,~, iter3, resvec3] = pcg(A, b, tol, maxit);
toc

L = speye(size(L));
tic
[~, resvec4, iter4] = mypcg(A, b, tol, maxit, L);
toc

figure(2)
semilogy(0:iter3, resvec3, 'g-o', 0:iter4, resvec4, 'r-*')
legend('Matlab CG', 'My PCG with L=Identity');
xlabel('Iterations');
ylabel('Residual Norm');

