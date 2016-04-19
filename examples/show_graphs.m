%
% Function to display the graph and the various clusters
% rw = 1 iff we want the vertices to be randomly perturbed
% to avoid collinear nodes
% if sw2 = 1, then show edges in cluster
%

function show_graphs(W,K,Xc,numblocks,nx,ny,sw2,rw)
  if rw == 1
      vertlist = drawgraph_v2(nx,ny); 
  else
      vertlist = drawgraph(W,nx,ny);
  end
  show_the_graph(W,nx,ny,vertlist)
  quit2 = 0;
     while quit2 ~= 1
      prompt = 'Cluster number (0 to end):';
      j = input(prompt);
      if j == 0 || j < 0 || j > K
         quit2 = 1;
      else
         fprintf('Number of nodes in cluster = %d \n',numblocks(j))
         [Wi, vx, vy] = drawblock_v2(W,Xc,vertlist,j);    
         [Wp,Wn] = splitW(Wi);
         figure
         hold on
         if sw2 == 1
            gplot(Wp,vertlist,'-b')
            gplot(Wn,vertlist,'-r')
         end
         plot(vertlist(:,1)', vertlist(:,2)','*r')
         plot(vx,vy,'*b')
         axis([0.7 nx + 0.3 0.7 ny + 0.3])
         hold off
      end
     end
end

