%
% Given an n x k matrix Z, this function tries to rescale the columns of
% Z to make the rows of Z  of unit length as much as possible
% If some entry in Y is negative, fail
% If all entries in Y are positive but one entry is very small, fail
%

function [Z2, Y, lam, flag, D] = normcol(Z)
scale = 1; flag = 0; tol = 10^(-14); 
n = size(Z,1); k = size(Z,2); 
lam = zeros(1,k); Z2 = zeros(n,k); D = zeros(1,k);
A = Z.*Z;
B = pinv(A);
% B
Y = B*scale*ones(n,1);
% Y
[minY,~] = min(Y);
if minY < 0
   flag = 1; return
else
   if minY < tol
      flag = 1; return
   else
     lam = Y.^(1/2); 
     Lam = diag(lam);
     Z2 = Z*Lam;
     NN = Z2*Z2';
     for i = 1:n
         D(i) = NN(i,i)^(1/2); 
     end
   end
end
end