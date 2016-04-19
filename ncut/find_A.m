%
% Given X and Z, looks for an invertible matrix A such that
% ||X - Z*A||_F is minimal
% Equivalently, -2tr(ZT*X*A) + tr(Z^T*Z*A^2) is minimized
% NNOTE -2tr(Z^T*X*A^T) + tr(Z^T*Z*A*A^T)
% Uses the explicit pseudo-inverse (Z'*Z)^{-1}*Z' 

function [A, fail] = find_A(X,Z,show1,show2)
   fail = 0; tol = 10^(15);
   k = size(Z,2);
   Y = Z'*X;
   % dy = det(Y);
   CY = cond(Y);
   ZZ = Z'*Z;
   A = inv(ZZ)*Y;
   ZA = Z*A;
   % d = det(A);
   C = cond(A);
   [U,S,V] = svd(A)
   sigk = S(k,k);
   if show2 == 1       
     % Z
     % X
     % ZA
     % fprintf('det(Y) = %d \n',dy)
      fprintf('Cond(Y)  = %d \n',CY)
     % fprintf('det(A) = %d \n',d)
      fprintf('Cond(A)  = %d \n',C)
      fprintf('Smallest singular value of A  = %d \n',sigk)
      fprintf('ZZ = Z^T*Z, Y = Z^T*X, and A in find_A \n')
    % ZZ
      Y
      A
   end
   if C > tol
       fail = 1;
     %  A = eye(k);
   end
   Diff1 = X - ZA;
      % Diff1
   NN1 = sqrt(trace(Diff1'*Diff1));
   if show2 == 1   
      fprintf('NN1 in find_A = %d \n',NN1)
   end
   Diff2 = X - Z;
   % Diff2
   NN2 = sqrt(trace(Diff2'*Diff2));
   if show2 == 1   
      fprintf('NN2 in find_A = %d \n',NN2)
   end 
end

