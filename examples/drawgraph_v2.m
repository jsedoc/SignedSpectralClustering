%
% To plot a graph whose nodes are the points on a nx x ny grid
% and the edges are specifided by an (nx x ny) x (nx x ny) adjacency matrix A
% This version perturbs the coordinates of the grid points by some small random
% vectors to avoid edge overlaps
%
%


function  vertlist = drawgraph_v2(nx, ny)
vertlist = zeros(nx*ny,2); scale = 0.5;
for j = 1:ny
     for i = 1:nx
         disturb = scale*(rand(1,2) - 0.5);
         vertlist((j-1)*nx + i,1) = i + disturb(1);
         vertlist((j-1)*nx + i,2) = j + disturb(2);
     end
end
end

