%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 7
% 
load("mat13041.rig");
A = spconvert(mat13041);
n = size(A, 1);
b1 = 1./sqrt(1:n);
x_exact = b1';
b = A * x_exact;
tol = 1e-12;
maxit = 550;
x0 = zeros(n,1);

setup.type = 'crout';
setup.droptol = 0.01;
[L,U] = ilu(A,setup);


k = 1;
for restart=[10,20,30,50]
    disp(restart)
    tic
    [x1,flag1,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,L,U);
    toc
    totalit = (iter(1)-1)*restart + iter(2)
    residuef = relres
    k = k+1;
    semilogy(0:totalit, resvec(1:totalit+1))
    hold on
    pause(1)
end


