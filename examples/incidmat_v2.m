%
% Given a symmetric matrix W with nonnegative entries
% and zero diagonal entries, computes a generalized
% incidence matrix for the oriented version of the underlying 
% graph of W obtained by directing edges from node i to node j iff i < j
% This version accomodates negative weights
% 

function B = incidmat_v2(W)
[nedges, edgelist] = elist(W);
nvert = size(W,1);
B = zeros(nvert,nedges);
for k = 1:nedges
    i = edgelist(k,1); j = edgelist(k,2);
    if W(i,j) >= 0
       B(i,k) = sqrt(W(i,j));
       B(j,k) = - B(i,k);
    else
       B(i,k) = sqrt(-W(i,j));
       B(j,k) = B(i,k);
    end
end
end