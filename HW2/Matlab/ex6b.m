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

iters = zeros(4,1);
residuef = zeros(4,1);
times = zeros(4,1);
labels = {'Matlab GMRES IC(0)';'My PCG IC(0)';'Matlab GMRES w/o Prec.';'My PCG with L=Identity'};

%% Original Exercise (only 1 iteration in both methods)

% Matlab GMRES Without Restart
tic
[~,~,~, iter1, resvec1] = gmres( A, b, restart, tol, maxit, L, L');
times(1) = toc;
totalit1 = (iter1(1)-1)*restart + iter1(2);
iters(1) = totalit1;
residuef(1) = resvec1(end);

% My PCG Implementation
tic
[~, resvec2, iter2] = mypcg(A, b, tol, maxit, L);
times(2) = toc;
iters(2) = iter2;
residuef(2) = resvec2(end);

% Residual Norm Plot
semilogy(0:totalit1, resvec1, 'g-o', 0:iter2, resvec2, 'r-*')
legend(labels(1:2));
xlabel('Iterations');
ylabel('Residual Norm');

% Results table
format shortEng
tab0 = [iters residuef times];
tab1 = tab0(1:2,:)
format default


%% Without Preconditioning (larger number of iterations)

% Matlab GMRES Without Restart
tic
[~,~,~, iter3, resvec3] = gmres(A, b, restart, tol, maxit);
times(3) = toc;
totalit3 = (iter3(1)-1)*restart + iter3(2);
iters(3) = totalit3;
residuef(3) = resvec3(end);

% My PCG Implementation
L = speye(size(L));
tic
[~, resvec4, iter4] = mypcg(A, b, tol, maxit, L);
times(4) = toc;
iters(4) = iter4;
residuef(4) = resvec4(end);

% Residual Norm Plot
figure(2)
semilogy(0:totalit3, resvec3, 'g-o', 0:iter4, resvec4, 'r-*')
legend(labels(2:3));
xlabel('Iterations');
ylabel('Residual Norm');

% Results table
format shortEng
tab2 = [iters residuef times]
format default


