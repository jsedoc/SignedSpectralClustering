function [U, D, V] = ransvd(A, k, p)
  if nargin < 2
      error(message('randsvd:TooFewInputs'));
  end    
  
  n = size(A,1);
  if k>n
      k=n;
  end
  if nargin < 3
      p = min(n-k,5);
  end
  
  Omega = randn(n, k+p);
  Y = A * Omega;
  [Q, ~] = qr(Y,0);
  B = Q' * A;
  [Uhat, D, V] = svd(B, 'econ');
  U = Q * Uhat;
end

