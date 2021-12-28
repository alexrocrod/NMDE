%% Exercise # 3 -  Numerical Solution of the Poisson Problem
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Input

% mesh
mesh = '0';
file = ['Input files\mesh' mesh '\mesh' mesh]; % 
coord = load([file '.coord']);
topol = load([file '.topol']);

% boundary conditions
% u = 0 for bound nodes
bound = load([file '.bound']);

% parameters
Ne = length(topol);
Nn = length(coord);
Rmax = 1e15;


%% Patter Creation for the Stifness Matrix
row1 = 1:length(topol);
row2 = [row1' row1' row1'];
row = reshape(row2',1,[]);
col = reshape(topol',1,[]);
A = sparse(row,col,1);
H = A' * A;

clear row1 row2 row col A;

%% Local Stifness Matrices
delta = ones(Ne,1);

for k = 1:Ne
    x = coord(topol(k,:),1);
    y = coord(topol(k,:),2);
    
    delta(k) = 0.5 * det([ones(3,1) coord(topol(k,:),:)]);
     
    valid_ijm = [1 2 3 1 2];
    b = zeros(1,3);
    c = zeros(1,3);
    for i1=1:3
        j1 = valid_ijm(i1+1);
        m1 = valid_ijm(i1+2);
%         a(i1) = x(j1)*y(m1) - x(m1)*y(j1);
        b(i1) = y(j1) - y(m1);
        c(i1) = x(m1) - x(j1);
    end

    B = b' * b;
    C = c' * c;
    Hloc = 0.25/delta(k) * (B + C); 

    % Algorithm 3.3
    for i = 1:3
        r1 = topol(k,i);
        for j = 1:3
            c1 = topol(k,j);
            H(r1,c1) = H(r1,c1) + Hloc(i,j);
        end
    end
end

%% Right-hand Side

f = zeros(Nn,1);
for i = 1:Nn
    x = coord(i,1);
    y = coord(i,2);
    tris = mod(find(topol==i),40);
    tris(tris==0) = 40; % fix changing 40, 80, 120 to 0
    f(i) = (-4 + 2*x*x + 2*y*y)*sum(delta(tris))/3;  
end

%% Boundary Conditions Enforcement

for i = bound(:,1)
    H(i,i) = Rmax;
    f(i) = Rmax*f(i);
end

%% Linear System Solution
tol = 1e-8;
maxit = 1000;

% Jacobi
M = sparse(diag(diag(H)));
tic
[x1, flag1, relres1, iter1, resvec1] = pcg(H, f, tol, maxit, M);
toc

% Choledsky
% L = ichol(H);
% tic
% [x2, flag2, relres, iter2, resvec2] = pcg(H, f, tol, maxit, L, L');
% toc



%% Error Computation

error = 0;
for i=1:Nn
    xsq = coord(i,1)^2;
    ysq = coord(i,2)^2;
    u_exact = xsq + ysq - xsq*ysq -1;
    u_exp = x1(i);

    tris = mod(find(topol==i),40);
    tris(tris==0) = 40; % fix changing 40, 80, 120 to 0
    suma = sum(delta(tris))/3; 

    res = (u_exp - u_exact)^2 * suma;
    error = error + res;
end
epsilon = sqrt(error)



