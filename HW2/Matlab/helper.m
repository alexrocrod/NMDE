% Extra helper script

%% Compute general data about a sparse matrix representing the finite difference Laplacian
close all
clear all

nxs = [102,202,402,802];
n = length(nxs);
conds = zeros(n,1);
sizes = zeros(n,1);
for i =1:n
    A = delsq(numgrid('S',nxs(i)));
    sizes(i) = size(A,1)
    conds(i) = condest(A)
end
sqc = sqrt(conds)
hs = 1./sizes;

%% Compute Error for ex3b
close all
clear all

bnorm = 811;
for o = 4:6
    o
    1/pk(o)/bnorm
end

function p = pk(k)
    lambs = [1,200,400,600,800,1000];
    p = max(chebyshevT(k-4,lambs));
end
