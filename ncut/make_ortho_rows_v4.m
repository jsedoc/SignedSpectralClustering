% Program to construct a k x k matrix R whose columns are rows
% of an n x k matrix A that are as orthogonal as possible
% Deletes zero rows, and normalizes the columns
%


function [R,fail] = make_ortho_rows_v4(A)
 tol = 10^(-15);
 fail = 0;
 %  deletes zero rows
 Izeros = sum(abs(A'));    % sum of the absolute values of the entries in each row 
 A = A(Izeros>0,:);
 n = size(A,1);
 k = size(A,2);
 if n < k
    fail = 1;
    R = eye(k); 
    return 
 end
 
 %R(:,1) = A(1,:)';
 % A2 = zeros(n-1,k);
 % A2 = A(2:n,:);
 
 % picks the first row
 index = 1;
 % pick the row of mid index
 % index = round(n/2);   
 R(:,1) = A(index,:)';
 A2 = zeros(n-1,k);
 A2(1:index-1,:) = A(1:index-1,:);
 A2(index:n-1,:) = A(index+1:n,:); 
 c = zeros(n-1,1);
 % A2
 for i = 2:k 
    c = c + abs(A2 * R(:,i - 1));
    [~, I] = min(c);
  %  I
  %  m
  %  cc= c';
  %  cc
  %  A2'
    R(:, i) = A2(I,:)';
    A3 = A2([1:I-1 I+1:n+1-i],:);
    %A3 = zeros(n-i,k);
    %A3(1:I-1,:) = A2(1:I-1,:);
    %A3(I:n-i,:) = A2(I+1:n+1-i,:);
    c2 = c([1:I-1 I+1:n+1-i],1);
    %c2 = zeros(n-i,1);
    %c2(1:I-1,1) = c(1:I-1,1);
    %c2(I:n-i,1) = c(I+1:n+1-i,1);
    %c = zeros(n-i,1);
    c = c2;
    %A2 = zeros(n-i,k);
    A2 = A3;
  %  A2'
 end
 NRc = zeros(k,k);
 RR = R'*R;
 for i = 1:k
     NRc(i,i) = RR(i,i)^(-1/2);
 end
 R = R*NRc';
end
