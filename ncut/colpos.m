%
% Multiply a column by -1 if necessary to make 
% the sum of entries nonnegative
%

function [Z2, RR] = colpos(Z)
N = size(Z,1); K = size(Z,2);
RR = eye(K);
colsum = ones(1,N)*Z;
Z2 = Z;
for i = 1:K
    if colsum(i) < 0
       Z2(:,i) = -Z(:,i);
       RR(i,i) = -1;        % records the sign flips in a diagonal matrix
    end
end
end
