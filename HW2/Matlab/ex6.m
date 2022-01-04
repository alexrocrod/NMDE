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
tol = 1e-8; %1e-10;
maxit = 550;
x0 = zeros(n,1);

setup.type = 'crout';
setup.droptol = 0.1;
[L,U] = ilu(A,setup);

[x, iter, resvec, flag] = myprecgmres(A, b, tol, maxit, x0, L, U);

norm(b-A*x)


%% c) from 3b)

clear all
close all

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
[x, flag, relres, iter1, resvec1] = gmres( A, b, restart, tol, maxit,L,L');
% [x, flag, relres, iter1, resvec1] = gmres( A, b, restart, tol, maxit)
toc
totalit = (iter1(1)-1)*restart + iter1(2);

% L = eye(size(L));
tic
[x, resvec, iter] = mypcg(A, b, tol, maxit, L);
toc

semilogy(0:totalit, resvec1, 'r-*',0:iter, resvec, 'g-+')
legend('Matlab GMRES' , 'My PCG');
xlabel('Iterations');
ylabel('Residual Norm');

