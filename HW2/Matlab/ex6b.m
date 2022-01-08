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
restart = n;

%% Original Exercise (only 1 iteration in both methods)

% Matlab GMRES Without Restart
tic
[~,~,~, iter1, resvec1] = gmres( A, b, restart, tol, maxit, L, L');
toc
totalit1 = (iter1(1)-1)*restart + iter1(2);

% My PCG Implementation
tic
[~, resvec2, iter2] = mypcg(A, b, tol, maxit, L);
toc

% Matlab PCG
tic
[~,~,~,  iter3, resvec3] = pcg(A, b, tol, maxit, L, L');
toc

% Residual Norm Plot
semilogy(0:totalit1, resvec1, 'g-o', 0:iter2, resvec2, 'r-*', 0:iter3, resvec3, 'b-+')
legend('Matlab GMRES IC(0)', 'My PCG IC(0)', 'Matlab PCG IC(0)');
xlabel('Iterations');
ylabel('Residual Norm');


%% Without Preconditioning (larger number of iterations)

% Matlab GMRES Without Restart
tic
[~,~,~, iter4, resvec4] = gmres(A, b, restart, tol, maxit);
toc
totalit4 = (iter4(1)-1)*restart + iter4(2);

% My PCG Implementation
L = speye(size(L));
tic
[~, resvec5, iter5] = mypcg(A, b, tol, maxit, L);
toc

% Matlab PCG
tic
[~,~,~,  iter6, resvec6] = pcg(A, b, tol, maxit);
toc

% Residual Norm Plot
figure(2)
semilogy(0:totalit4, resvec4, 'g-o', 0:iter5, resvec5, 'r-*', 0:iter6, resvec6, 'b-+')
legend('Matlab GMRES w/o Prec.', 'My PCG with L=Identity', 'Matlab CG');
xlabel('Iterations');
ylabel('Residual Norm');


