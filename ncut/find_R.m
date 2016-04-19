%
%  Function to find a matrix R in O(K) that minimizes
%  ||X - ZR||_F, given X and Z 
%

function R = find_R(X,Z)
  A = Z'*X;
  [U, ~, V] = svd(A);
  R = U*V';
end