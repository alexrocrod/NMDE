%% Exercise # 3 -  Numerical Solution of the Poisson Problem
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Input
prec = 'C'; % J -> use Jacobi, C -> use Cholesky

% PCG Parameters
tol = 1e-8;
maxit = 200;

% Folder to save results
savefolder = ['Results\' prec '\'];
mkdir(savefolder(1:end-1))

resvecs = zeros(maxit,5);
iters = zeros(1,5);
Tsols = zeros(1,5);
epsilons = zeros(1,5);

for imesh=0:4
    % Load mesh files
    mesh = num2str(imesh);
    file = ['Input files\mesh' mesh '\mesh' mesh]; 
    coord = load([file '.coord']);
    topol = load([file '.topol']);
    bound = load([file '.bound']);
    
    % parameters
    Ne = length(topol);
    Nn = length(coord);
    Rmax = 1e15;
    
    %% Patter Creation for the Stifness Matrix
    row1 = 1:Ne; % range 
    row2 = [row1' row1' row1']; % 3 ranges as colloums
    row = reshape(row2',[],1); % row2 as a column vector
    col = reshape(topol',[],1); % topol as a column vector
    A = sparse(row,col,1); 
    H = A' * A;
    H = H * 0; % clean coefficients form H
    
    clear row1 row2 row col A;
    
    %% Stifness Matrix
    [H, delta] = computeStiff(H, topol, coord);
    
    %% Right-hand Side
    
    f = zeros(Nn,1);
    for i = 1:Nn
        x = coord(i,1);
        y = coord(i,2);
    
        els = mod(find(topol==i),Ne); % elements that have that node
        els(els==0) = Ne; % fix changing k*Ne to 0
    
        f(i) = (-4 + 2*x^2 + 2*y^2) * sum(delta(els)) / 3;  
    end
    
    %% Boundary Conditions Enforcement
    
    for i = bound(:,1)'
        H(i,i) = Rmax; % substitute
    end
    
    %% Linear System Solution
    
    if prec == 'J' % Jacobi
        M = sparse(diag(diag(H)));
        tic
        [u, flag, relres, iter, resvec] = pcg(H, f, tol, maxit, M);
        Tsol = toc;
    elseif prec == 'C' % Cholesky
        L = ichol(H);
        tic
        [u, flag, relres, iter, resvec] = pcg(H, f, tol, maxit, L, L');
        Tsol = toc;
    else
        disp('Invalid value for the preconditioner id., valid values: J and C.')
        return;
    end

    %save results
    resvecs(1:iter+1,imesh+1) = resvec;
    iters(imesh+1) = iter;
    Tsols(imesh+1) = Tsol;
    
    %% Error Computation
    
    error = 0;
    for i=1:Nn
        % Compute analytical solution
        xi = coord(i,1);
        yi = coord(i,2);
        u_exact = xi^2 + yi^2 - xi^2 * yi^2 - 1; % analytical solution
        
        % Compute surface measure
        els = mod(find(topol==i),Ne);
        els(els==0) = Ne; % fix bad change of k*Ne to 0
    
        res = (u(i)- u_exact)^2 * sum(delta(els))/3;
        error = error + res;
    end
    epsilons(imesh+1) = sqrt(error);
    
    
end
%% Save Results


% Convergence Plot
semilogy(0:iters(1),resvecs(1:iters(1)+1,1),'r-*', ...
    0:iters(2),resvecs(1:iters(2)+1,2),'g-o', ...
    0:iters(3),resvecs(1:iters(3)+1,3),'b-+', ...
    0:iters(4),resvecs(1:iters(4)+1,4),'k-x', ...
    0:iters(5),resvecs(1:iters(5)+1,5),'c-h')
legend('Level 0', 'Level 1', 'Level 2', 'Level 3', 'Level 4')
axis([0 1.1*max(iters) 0 max(max(resvecs))])
xlabel('Iterations');
ylabel('Residual Norm');

% Save figure
f = gcf;
exportgraphics(f, [savefolder 'Convergence.png'])

% Save epsilon and computational time in a .txt file
save([savefolder 'res.txt'],'epsilons','Tsols','-ascii');

%% Functions
% Compute Stifness Matrix
function [H, delta] = computeStiff(H, topol, coord)
    
    Ne = length(topol);
    delta = ones(Ne,1);
    
    for k = 1:Ne
        x = coord(topol(k,:),1);
        y = coord(topol(k,:),2);
        
        % Surface Meause of each element
        delta(k) = 0.5 * det([ones(3,1) coord(topol(k,:),:)]);
         
        valid_ijm = [1 2 3 1 2]; % simple rotation of indices
        % compute b and c
        b = zeros(1,3);
        c = zeros(1,3);
        for i1=1:3
            j1 = valid_ijm(i1+1);
            m1 = valid_ijm(i1+2);
            b(i1) = y(j1) - y(m1);
            c(i1) = x(m1) - x(j1);
        end
        % Hloc Computation
        B = b' * b; % [bibi, bibj, ...]
        C = c' * c;
        Hloc = 1/(4*delta(k)) * (B + C); 
    
        % Algorithm 3.3 to assemble H
        for i = 1:3
            r1 = topol(k,i);
            for j = 1:3
                c1 = topol(k,j);
                H(r1,c1) = H(r1,c1) + Hloc(i,j);
            end
        end
    end
end
