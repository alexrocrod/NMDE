
% close all
% clear all

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
