% code from 
% https://www.cs.purdue.edu/homes/dgleich/demos/matlab/spectral/spectral.html

clear A;

n = 10;
gs2 = floor(0.4*n);
gs3 = floor(0.2*n);

%x = randperm(n);
x = 1:n;

group1 = x(1:gs2);
group2 = x(gs2 + 1:gs2 + gs3);
group3 = x(gs2 + gs3 + 1:end);

p_group1 = 0.8;
p_group2 = 0.7;
p_group3 = 0.7;
p_between = 0.1;

A(group1, group1) = rand(gs2,gs2)     < p_group1;
A(group2, group2) = rand(gs2 - gs3,gs2 - gs3) < p_group2;
A(group3, group3) = rand(n - (gs2 + gs3) ,n - (gs2 + gs3)) < p_group3;
A(group1, group2) = rand(size(group1,2), size(group2,2))  < p_between;
A(group1, group3) = rand(size(group1,2), size(group3,2))  < p_between;
A(group2, group3) = rand(size(group2,2), size(group3,2))  < p_between;


% remove isolates nodes.
D = sum(abs(A));
A = A(D>0,D>0);

A = triu(A,1);
A = A + A';
spy(A);
A