%
% Tries to make the rows of Z2 as close to rowsum as possible
% Using the least squares method. It fails if some
% lambda = 0 in lam
%

function [Zr, Lam, fail] = rowsumto1(Z2,rowsum,show2)
 tol = 10^(-10);
 fail = 0;
 K = size(Z2,2); m = size(Z2,1);
 Z2Z2 = Z2'*Z2; DZ2inv = zeros(1,K);
 for i = 1:K
     DZ2inv(i)= Z2Z2(i,i)^(-1);
 end
 Z2Z2inv = diag(DZ2inv); 
 lam = rowsum*Z2Z2inv*Z2'*ones(m,1);
 if show2 == 1
    fprintf('lam in rowtosum1  \n')
    lam
 end
 [m,~] = min(abs(lam));
 if m < tol
    fail = 1; 
 end
    Lam = diag(lam');
    Zr = Z2*Lam;
end

