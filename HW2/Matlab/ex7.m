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

restarts = [10,20,30,50];
n = length(restarts);
resvecs = zeros(n,maxit*50);
totalits = zeros(n,1);
residuef = zeros(n,1);
times = zeros(n,1);

idx = 1;
for restart=restarts
    disp(restart)
    tic
    [x1,flag1,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,L,U);
    times(idx) = toc;
    totalit = (iter(1)-1)*restart + iter(2);
    totalits(idx) = totalit;
    residuef(idx) = relres;
    resvecs(idx,1:totalit+1) = resvec(1:totalit+1)';

    semilogy(0:totalit, resvec(1:totalit+1))
    hold on
    pause(1)
    idx = idx+1;
end

figure(2)
for i=1:n
    semilogy(0:totalits(i), resvecs(i,1:totalits(i)+1))
    
    hold on
end
axis([0 max(totalits) min(min(resvecs~=0)) max(max(resvecs))])
legend('Restart=10','Restart=20','Restart=30','Restart=50')

format shortEng
tab = [restarts' totalits residuef times]
format default

