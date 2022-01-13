%% Exercise # 3 -  Numerical Solution of the Poisson Problem
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Input
prec = 'J'; % J -> use Jacobi, C -> use Choledsky
mesh = '0'; % mesh level of refinement

% Folder to save results
savefolder = ['Results\' prec  '_' mesh '\'];
mkdir(savefolder(1:end-1))

% Load from mesh files
file = ['Input files\mesh' mesh '\mesh' mesh]; 
coord = load([file '.coord']);
topol = load([file '.topol']);

% boundary conditions
% u = 0 for boundary nodes
bound = load([file '.bound']);

% parameters
Ne = length(topol);
Nn = length(coord);
Rmax = 1e15;

%% Patter Creation for the Stifness Matrix
row1 = 1:length(topol); % range 
row2 = [row1' row1' row1']; % 3 ranges as colloums
row = reshape(row2',1,[]); % row2 as a row vector
col = reshape(topol',1,[]); % topol as a row vector
A = sparse(row,col,1); 
H = A' * A;

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

for i = bound(:,1)
%     H(i,:) = 0; destroys simetry
    H(i,i) = Rmax*H(i,i); % works
%     H(i,i) = Rmax; % substitute -> gives ichol(H) error
%     f(i) = Rmax*f(i); % not needed in our case
end

%% Linear System Solution
% PCG Parameters
tol = 1e-8;
maxit = 1000;

if prec == 'J' % Jacobi
    M = sparse(diag(diag(H)));
    tic
    [x, flag, relres, iter, resvec] = pcg(H, f, tol, maxit, M);
    Tsol = toc
elseif prec == 'C' % Choledsky
    L = ichol(H);
    tic
    [x, flag, relres, iter, resvec] = pcg(H, f, tol, maxit, L, L');
    Tsol = toc
else
    disp('Invalid value for the preconditioner id, valid values: J and C.')
end

% Convergence Plot
semilogy(0:iter,resvec,'r-*')
xlabel('Iterations');
ylabel('Residual Norm');

% Save figure
f = gcf;
exportgraphics(f, [savefolder 'Convergence.png'])


%% Error Computation

error = 0;
for i=1:Nn
    % Compute analytical solution
    xsq = coord(i,1)^2;
    ysq = coord(i,2)^2;
    u_exact = xsq + ysq - xsq*ysq - 1; % analytical solution
    
    % Compute surface measure
    els = mod(find(topol==i),Ne);
    els(els==0) = Ne; % fix bad change of k*Ne to 0, k=1,2,...
    surf = sum(delta(els))/3; % surface measure

    res = (x(i)- u_exact)^2 * surf;
    error = error + res;
end
epsilon = sqrt(error)

% save epsilon and computational time in a .txt file
save([savefolder 'res.txt'],'epsilon','Tsol','-ascii');


%% Functions
% Compute Stifness Matrix
function [H, delta]= computeStiff(H, topol, coord)
    
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
