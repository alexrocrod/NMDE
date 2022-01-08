%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 7
A = load("mat13041.rig");
A = spconvert(A);
n = size(A, 1);
b1 = 1./sqrt(1:n);
x_exact = b1';
b = A * x_exact;
tol = 1e-12;
maxit = 550;
x0 = zeros(n, 1);

% Preconditioning
setup.type = 'crout';
setup.droptol = 0.01;
[L,U] = ilu(A, setup);

restarts = [10, 20, 30, 50];
n = length(restarts);
resvecs = zeros(n, maxit * 50);
totalits = zeros(n, 1);
residuef = zeros(n, 1);
times = zeros(n, 1);

idx = 1;
for restart=restarts
    % run method
    tic
    [~, ~, relres, iter, resvec] = gmres(A, b, restart, tol, maxit, L, U);
    times(idx) = toc;

    % save results
    totalit = (iter(1) - 1) * restart + iter(2);
    totalits(idx) = totalit;
    residuef(idx) = resvec(end);
    resvecs(idx, 1:totalit+1) = resvec(1:totalit+1)';
    
    % display results
    fprintf('Restart = %d -> Iterations: %d, Relative Residue: %d, Time: %d s\n', ...
        restart, totalit, residuef(idx), times(idx))

    idx = idx + 1;
end

% Residual Norm vs Iterations Plot
linestyles = {'r-*', 'g-o', 'b-+', 'k-h'};
for i=1:n
    semilogy(0:totalits(i), resvecs(i, 1:totalits(i)+1), linestyles{i})
    hold on
end
axis([0 max(totalits) 0 max(max(resvecs))])
legend('Restart=10','Restart=20','Restart=30','Restart=50')

% Results table
format shortEng
tab = [restarts' totalits residuef times]
format default

