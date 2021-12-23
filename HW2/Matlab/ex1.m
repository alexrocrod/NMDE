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

% Matlab PCG without preconditioning
tic
[x, flag, relres, iter1, resvec1] = pcg( A, b, tol, maxit); 
toc
% Matlab PCG
tic
[x, flag, relres, iter2, resvec2] = pcg( A, b, tol, maxit, L, L');
toc

% my implementation
tic
[x, resvec3, iter3] = mypcg(A, b, tol, maxit, L);
toc
semilogy(0:iter1, resvec1, 'r-*', 0:iter2, resvec2,'g-o', 0:iter3, resvec3,'b-+')


%% Function
function [x, resvec, iter] = mypcg(A, b, tol, maxit, L)
    M1 = L;
    M2 = L';
    M = M1 * M2;
    resvec = zeros(maxit+1, 1);
        
    x = zeros(length(b),1);
    r = b - A * x;
    z = M^(-1) * r;
    p = z;
    rsold = r' * z;

    
    for iter = 1:maxit 
        Ap = A * p;
        alpha = rsold / (Ap' * p);
        x = x + alpha * p;
        r = r - alpha * Ap;

        v = M1\r;
        z = M2\v;  % w = M2\v;

        rsnew = r' * z;

        resvec(iter) = sqrt(rsnew);
        if sqrt(rsnew) < tol
              break
        end

        p = r + (rsnew / rsold) * p;
        rsold = rsnew; 
    end
    resvec = resvec(1:iter+1,:);

end

