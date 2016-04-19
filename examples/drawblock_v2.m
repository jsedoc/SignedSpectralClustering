%
% Function to draw the induced graph whose set of nodes
% belongs to the jth block of the partition
% This version uses the list of coordinates of the nodes, vertlist
%

function [Wi, vx, vy] = drawblock_v2(W,Xc,vertlist,j)
N = size(Xc,1); K = size(Xc,2);
vx = zeros(1); vy = zeros(1);
j1 = 0;
for k = 1:N
    if Xc(k,j) > 0 
        j1 = j1 + 1; vx(j1) =  vertlist(k,1);
        vy(j1) = vertlist(k,2);
    end
end
Wi = W;
for i = 1:N
    if Xc(i,j) == 0
       Wi(i,:) = 0; Wi(:,i) = 0;
    end
end
end