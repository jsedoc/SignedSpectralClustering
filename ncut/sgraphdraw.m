%
% Function to draw a graph or signed graph using the signed Laplacian
% This version also draws the graph corresponding to the matrix of 
% absolute values of the weight matrix 
%

function L = sgraphdraw(W, names)
m = size(W,1);   tol = 10^(-15); % tolerance to decide when an eigenvalue is 0
minval = min(min(W));   % finds smallest entry in W
if minval < 0
     fprintf('The graph is a signed graph  \n')
else
     fprintf('The graph is an unsigned graph  \n')
end
D = diag(sum(abs(W)));   %  computes degree matrix of absolute values
figure(1)
%colormap(gray)
imagesc(50*W)
L = D - W;             % standard signed Laplacian
% L
Wa = abs(W);           % matrix of absolute values of the entries in W
La = D - Wa;           % Laplacian of Wa
prompt = 'Method (1 = svd(L), 2 = eig(L)):';
mm = input(prompt);
  fprintf('Number of nodes = %d \n',m)
  if mm == 1
     [U,S,V] = svd(L); % SVD of Laplacian
     [Ua,Sa,Va] = svd(La); % SVD of Laplacian of La 
     vertlista = Ua(:,[m-1 m-2]);
     S 
     U
     lamb1 = S(m,m);
  else
     [v, e] = eig(L); 
     [va,ea] = eig(La);
     vertlista = va(:,[2 3]);
      e
      v
     lamb1 = e(1,1);
  end
 
  sw = 1;
  fprintf('smallest eigenvalue of L = %d \n',lamb1)
  if minval < 0
     if lamb1 < tol 
        fprintf('The graph is balanced  \n')
        prompt = 'Balanced graph plotted as a bipartite graph (no = 0):';
        sw = input(prompt);
     else
        fprintf('The graph is unbalanced  \n')
     end
  end
  
  if minval < 0       % the graph has some negative edge
      if sw == 1      % the graph is unbalanced or plotted as bipartite
         if mm == 1
            vertlist = U(:,[m m-1]);
         else
            vertlist = v(:,[1 2]);  
         end
      else            % the graph is balanced and not treated as bipartite 
         if mm == 1
            vertlist = U(:,[m-1 m-2]);
         else
            vertlist = v(:,[2 3]);  
         end
      end
  else                % the graph has no negative edge
      if mm == 1
        vertlist = U(:,[m-1 m-2]);
      else
        vertlist = v(:,[2 3]);  
      end
  end
  [Wp,Wn] = splitW(W);   %  splits W into its positive part and its negative part
  if minval < 0
     figure
     hold on
     gplot(Wp,vertlista,'-b')
     gplot(Wn,vertlista,'-r')
     gplot(W,vertlista,'*r')
     A = vertlista
     if nargin > 1
      for li=1:size(A,1) 
        text( A(li,1), A(li,2), sprintf('    %s', char(names(li))) );
      end
      hold off
     end
  end
  figure
  hold on
  gplot(Wp,vertlist,'-b')
  gplot(Wn,vertlist,'-r')
  gplot(W,vertlist,'*r')
  A = vertlist
  if nargin > 1
    for li=1:size(A,1)  
      text( A(li,1), A(li,2), sprintf('    %s', char(names(li))) );
    end
  end
  hold off
end

