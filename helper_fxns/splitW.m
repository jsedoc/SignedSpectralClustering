%
%  Takes a symmetric weight matrix W with positive and negative
%  weights and retunrs the matrices Wp and Wn corresponding to the
%  positive weights and negative weights (with positive weights)
%

function [Wp, Wn] = splitW(W)
m = size(W,1);
Wp = W;
for i = 1:m
    for j = 1:m
        if Wp(i,j) < 0
           Wp(i, j) = 0;
        end
    end
end
Wn = Wp - W;
end