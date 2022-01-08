%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 8 
A = load("ML_laplace.mtx");
A = spconvert(A);
nA = size(A, 1);
x_exact = ones(nA,1);
b = A * x_exact;
tol = 1e-8;
maxit = 550;
x0 = zeros(nA,1);

restart = 50;

tols = [2e-2,1e-2,3e-3,1e-3,1e-4,1e-5];
n = length(tols);
totalits = zeros(n,1);
residuef = zeros(n,1);
tsol = zeros(n,1);
tprec = zeros(n,1);
rho = zeros(n,1);
resvecs = zeros(n,maxit*restart);

idx = 1;
for dtol = tols    
    % Preconditioner Computation
    tic
    setup.type = 'crout';
    setup.droptol = dtol;
    [L,U] = ilu(A,setup);
    tprec(idx) = toc;
    
    % run method
    tic
    [x1,flag1,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,L,U);
    tsol(idx) = toc;

    % sabe results
    totalit = (iter(1)-1)*restart + iter(2);
    totalits(idx) = totalit;
    residuef(idx) = resvec(end);
    rho(idx) = (nnz(L) + nnz(U) - n)/nnz(A);
    resvecs(idx,1:totalit+1) = resvec(1:totalit+1)';

    % display results
    fprintf(['Droptol.: %d -> Iterations: %d, : Time Prec.: %d, Time Sol.: %d, ' ...
        'Residue: %d , Rho: %d \n'], dtol, totalit, tprec(idx), ...
        tsol(idx), residuef(idx), rho(idx))

    idx = idx + 1;
end

% Residual Norm Plot
linestyles = {'r-*', 'g-o', 'b-+', 'k-h', 'c-s', 'y-d'};
for i=1:n
    semilogy(0:totalits(i), resvecs(i,1:totalits(i)+1), linestyles{i})
    hold on
end
axis([0 1.1*max(totalits) 0 max(max(resvecs))])
legend('tol=2e-2','tol=1e-2','tol=3e-3','tol=1e-3','tol=1e-4','tol=1e-5')

% Results Table
ttotal = tprec + tsol;
format shortEng
tab = [tols' totalits tprec tsol ttotal residuef rho]
format default

