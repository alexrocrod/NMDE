%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 2
% 
nxs = [102 202 402 802];
ifig = 1;
for nx=nxs
    fprintf('nx = %d\n',nx);
    A = delsq(numgrid('S',nx));
    n = size(A, 1);
    b1 = 1./sqrt(1:n);
    x_exact = b1';
    b = A * x_exact;
    tol = 1e-8;
    maxit = 5000;
    
    % a)
    [x, flag, relres, iter1, resvec1] = pcg(A, b, tol, maxit); 


    % b)
    L = ichol(A);
    [x, flag, relres, iter2, resvec2] = pcg(A, b, tol, maxit, L, L');

    % c)
    opts.type = 'ict';
    opts.droptol = 1e-2;
    L = ichol(A,opts);
    [x, flag, relres, iter3, resvec3] = pcg(A, b, tol, maxit, L, L'); 

    % d)
    opts.type = 'ict';
    opts.droptol = 1e-3;
    L = ichol(A,opts);
    [x, flag, relres, iter4, resvec4] = pcg( A, b, tol, maxit, L, L');

    % error plots
    figure(ifig);
    semilogy(0:iter1, resvec1, 'r-*', 0:iter2, resvec2,'g-o', 0:iter3, resvec3,'b-+', 0:iter4, resvec4,'k-o')
    legend('No preconditioner' , 'IC(0)', 'droptol=1e-2','droptol=1e-3');
    xlabel('Iterations');
    ylabel('Residual Norm');
    ifig = ifig + 1;
    pause(1)

end

%Produce a Table with the values of h and the number of (P)CG iterations:
%one row for each nx value

% h vem de onde??

