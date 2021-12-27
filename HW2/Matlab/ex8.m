%% Exercise # 2 -  Iterative Methods For Linear Systems
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Question 8
% 
A1=load("ML_laplace.mtx");
A = spconvert(A1);
n = size(A, 1);
x_exact = ones(n,1);
b = A * x_exact;
tol = 1e-12;
maxit = 550;
x0 = zeros(n,1);

