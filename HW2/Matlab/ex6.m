%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 6
A = load("mat13041.rig");
A = spconvert(A);
n = size(A, 1);
b1 = 1./sqrt(1:n);
x_exact = b1';
b = A * x_exact;
tol = 1e-10;
maxit = 550;
x0 = zeros(n,1);


% Preconditioner
setup.type = 'crout';
setup.droptol = 0.01;
tic
[L,U] = ilu(A,setup);
toc

% Matlab GMRES Without Restart
restart = n;
tic
[x1, ~, ~, iter1, resvec1] = gmres(A, b, restart, tol, maxit,L,U);
toc
totalit = (iter1(1)-1)*restart + iter1(2);
trueres1 = norm(b-A*x1);
fprintf('Matlab GMRES -> Final Residue: %d, True Residue: %d\n', resvec1(end), trueres1)

% My GMRES Implementation
tic
[x2, iter2, resvec2, ~] = myprecgmres(A, b, tol, maxit, x0, L, U);
toc
trueres2 = norm(b-A*x2);
fprintf('My GMRES -> Final Residue: %d, True Residue: %d\n ', resvec2(end), trueres2)

% Residual Norm Plot
semilogy(0:totalit, resvec1, 'r-*', 0:iter2, resvec2, 'g-+')
legend('Matlab GMRES', 'My Prec. GMRES');
xlabel('Iterations');
ylabel('Residual Norm');



