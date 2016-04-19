function [ XXc ] = idx2xxc( IDX, num_clusters )
   if nargin <2
        num_clusters = max(IDX);
   end
   XXc = zeros(numel(IDX), num_clusters);
   for i =1:numel(IDX) 
       XXc(i,IDX(i)) = 1;
   end
end

