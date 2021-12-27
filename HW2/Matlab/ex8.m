%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 8
% 
A1=load("ML_laplace.mtx");
A = spconvert(A1);
n = size(A, 1);
x_exact = ones(n,1);
b = A * x_exact;
tol = 1e-12;
maxit = 550;
x0 = zeros(n,1);

restart = 50;

tols = [2e-2,1e-2,3e-3,1e-3,1e-4,1e-5];
for dtol = tols
    disp(dtol);

    tic
    setup.type = 'crout';
    setup.droptol = dtol;
    [L,U] = ilu(A,setup);
    toc
    
    tic
    [x1,flag1,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,L,U);
    toc

    totalit = (iter(1)-1)*restart + iter(2)
    residuef = relres
    rho = (nnz(L) + nnz(U) - n)/nnz(A)


    semilogy(0:totalit, resvec(1:totalit+1))
    hold on
    pause(0.01)

end

