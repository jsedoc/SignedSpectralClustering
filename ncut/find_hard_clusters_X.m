%
%
% Function to find a "discrete" partition matrix X that
% minimizes ||X - Z*R||_F, given Z and R in O(K) 
% We make sure that ||X||_F = ||Z||_F.
% This version also makes the sum of every column of Z*R positive 
% and picks the X with the smallest norm ||X - Z*R||_F
%

function [X,Rn] = find_hard_clusters_X(Z,R,a,show1,show2)
  ZR1 = Z*R;   %  an N x K matrix
  K = size(ZR1,2);
  Rn = eye(K);
  if show1 == 1
     fprintf('Z*R in find_hard_clusters_X \n')
     ZR1
  end
  [ZR2, Rn] = colpos(ZR1);   %  to make the sum of each column nonnegative
  if show2 == 1
    % fprintf('Z*R and Rn after colpos in find_hard_clusters_X \n')
     fprintf('Rn after colpos in find_hard_clusters_X \n')
    % ZR2
     Rn
  end
  for st = 1:2
      if st == 1
          ZR = ZR1;
      else 
          ZR = ZR2;
      end
      N = size(ZR,1); 
      X1 = zeros(N,K);   
      [~, J] = max(ZR');  % max entries in each row and their first index
       for i = 1:N
           X1(i,J(i)) = 1;
       end
       % X1
       X1sum = sum(X1);   % vector of the sums of the  columns
        % X1sum
      [~, rowindex] = max(X1);  % nonzero entries in each column and their first index
      [~, k] = max(X1sum);  % maximum number of 1's in each columns and first index 
                            % of such a column
      %  If some column j has all zeros, then reallocate a 1 from
      %  a column k with a maximum number of 1's; the first row with a nonzero
      %  in column k is rowindex(k)
      %
      for j = 1:K
          if X1sum(j) == 0
             if show1 ==1 
                fprintf('X1 in find_hard_clusters_X when there is a zero column, and st = %d \n',st)
             end
             X1(rowindex(k),k) = 0; 
             X1(rowindex(k),j) = 1;
             % X1
            X1sum = sum(X1);
            [~, rowindex] = max(X1);
            [~, k] = max(X1sum);
          end
      end    
      if show1 == 1
         fprintf('X1 in find_hard_clusters_X, and st = %d \n',st)
         X1
      end
      if st == 1
         XX1 = a*X1;   %  Normalizes X1 so that ||X||_F = ||Z||_F 
      else
         XX2 = a*X1;
      end
  end   
  N1 = sqrt(trace((XX1 - ZR1)'*(XX1 - ZR1)));
  N2 = sqrt(trace((XX2 - ZR2)'*(XX2 - ZR2)));
  if show1 == 1
     fprintf('||XX1 - Z*R||_F = %d \n',N1)
     fprintf('||XX2 - Z*R||_F = %d \n',N2)
  end
  if N1 > N2     % picks the X that yields the smallest norm ||X - Z*R||_F
        X = XX2;
  else
        X = XX1;
  end
  if show1 == 1
     X
  end
end