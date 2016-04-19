%
% To display the graph
%
function show_the_graph(W,nx,ny,vertlist) 
  [Wp,Wn] = splitW(W);   %  splits W into its positive part and its negative part
  figure
  hold on
  gplot(Wp,vertlist,'-*b')
  gplot(Wn,vertlist,'-*r')
  axis([0.7 nx + 0.3 0.7 ny + 0.3]) 
  hold off
end

