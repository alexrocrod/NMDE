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
lnx = length(nxs);
ifig = 1;
iters = zeros(lnx,lnx);

for nx=nxs
    fprintf('nx = %d\n',nx);
    A = delsq(numgrid('S',nx));
    n = size(A, 1);
    b1 = 1./sqrt(1:n);
    x_exact = b1';
    b = A * x_exact;
    tol = 1e-8;
    maxit = 5000;

    cond = condest(A) % condition number estimation
    
    % a) No Preconditioner
    tic
    [~, ~, ~, iter1, resvec1] = pcg(A, b, tol, maxit); 
    NoPrec = toc % computational time


    % b) PCG with IC(0)
    L = ichol(A);
    tic
    [~, ~, ~, iter2, resvec2] = pcg(A, b, tol, maxit, L, L');
    IC0 = toc % computational time
    

    % c) PCG with IC and droptol = 1e-2
    opts.type = 'ict';
    opts.droptol = 1e-2;
    L = ichol(A,opts);
    tic
    [~, ~, ~,  iter3, resvec3] = pcg(A, b, tol, maxit, L, L'); 
    IC2 = toc % computational time


    % d) PCG with IC and droptol = 1e-3
    opts.type = 'ict';
    opts.droptol = 1e-3;
    L = ichol(A,opts);
    tic
    [~, ~, ~,  iter4, resvec4] = pcg(A, b, tol, maxit, L, L');
    IC3 = toc % computational time

    % Residual Norm Plots
    figure(ifig);
    semilogy(0:iter1, resvec1, 'r-*', 0:iter2, resvec2, 'g-o', 0:iter3, ...
        resvec3, 'b-+', 0:iter4, resvec4, 'k-h')
    legend('No preconditioner', 'IC(0)', 'droptol=1e-2', 'droptol=1e-3');
    xlabel('Iterations');
    ylabel('Residual Norm');
    ifig = ifig + 1;
    if ifig <= lnx, pause(1), end
end


