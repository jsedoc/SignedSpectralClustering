%
% Create a list of directed edges  corresponding to the 
% positive entries in a symmetric matrix W
% with nonnegative entries (and zero diagonal entries)
%

function [c, edgelist] = elist(W)
   n = size(W,1);
   c = 0;
   edgelist = zeros(1,2);
   for i = 1:n
       for j = 1:n
           if i < j && W(i,j) ~= 0
              c = c + 1;
              edgelist(c,1) = i;
              edgelist(c,2) = j;
           end
       end
   end
end