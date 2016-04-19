function [ IDX, XXc ] = sncut(W, k, varargin)
%SNCUT Signed normalized cut method for clustering.
%   IDX = SNCUT(W, K) partitions the nodes described by the symmetric 
%         m-by-m degree matrix W into K clusters. 
%         SNCUT returns a m-by-1 vector IDX containing the cluster indices 
%         of each point.

%   SNCUT treats NaNs as missing data, and ignores any rows of X that
%   contain NaNs.

% References: 
%             
%           
%   [1] Gallier, J. H. Spectral Theory of Unsigned and Signed Graphs
%       Applications to Graph Clustering: a Survey. 
%       http://www.cis.upenn.edu/~jean/spectral-graph-notes.pdf


if nargin < 2
    error(message('sncut:TooFewInputs'));
end

wasnan = any(isnan(W),2);
hadNaNs = any(wasnan);
if hadNaNs
    error('sncut:W contains NaNs');
    %warning(message('sncut:MissingDataRemoved'));
    %W = W(~wasnan,:);
end

if max(max(abs(W-W'))) > 1e-10
    error('sncut:W not symmetric');
end

if max(diag(W)) > 1e-10
    error('sncut:W has non-zero diagonal entries');
end

% m vertices
m = size(W,1);

if k<2 || k >m
    error('sncut:Incorrect number of K cluster');
end

pnames = {   'maxiters'  'initZmethod' 'findAmethod' 'threshold' 'display' 'type'};
dflts =  {   50           1            1             1e-14        0        'hard'};
[maxiter,initZmethod,findAmethod,tol,display, type] ...
    = parseArgs(pnames, dflts, varargin{:});

if strcmp(type, 'hard')
    find2_X = @find_hard_clusters_X;
end
if strcmp(type, 'flex')
    find2_X = @find_flex_K_clusters_X;
end
if strcmp(type, 'soft')
    find2_X = @find_soft_clusters_X;
end



minval = min(min(W));   % finds smallest entry in W
if minval < 0
     fprintf('The graph is a signed graph  \n')
else
     fprintf('The graph is an unsigned graph  \n')
end

Dv = sum(abs(W));       % vector of degrees
if sum(Dv < tol)>0
    error('sncut:W has 0 column which means that the graph has isolated nodes');
end

D = diag(Dv);           %  computes degree matrix of absolute values
degv = D*ones(m,1); d = ones(m,1)'*degv; % vector of degrees and total volume
Dnhalf = diag(Dv.^(-1/2));



L = D - W;             % standard signed Laplacian
Ls = Dnhalf*L*Dnhalf;  % normalized  Laplacian

%%%% NOTE: for speed we are not using the generalized incidence matrix
%B = incidmat_v2(W);    % generalized incidence matrix of W
% B
%n = size(B,2);
%Bs = Dnhalf*B;         % normalized incidence matrix
Bs = 0;
%if issparse(Ls)
%    [U5,S5,V5] = svds(2*eye(m) - Ls,k);        % SVD of 2*I - Ls
%else
%    [U5,S5,V5] = svd(2*eye(m) - Ls);        % SVD of 2*I - Ls
%end
[U5,S5,V5] = fsvd(2*eye(m) - Ls,k);

lamb1 = S5(1,1);
  fprintf('lamb1 = %d \n',lamb1)
  fprintf('smallest eigenvalue of Ls = %d \n',2-lamb1)
  if minval < 0
     if abs(2-lamb1) < tol 
        fprintf('The graph is balanced  \n')
     else
        fprintf('The graph is unbalanced  \n')
     end
  end
  
  %  Prompts for method to compute eigenvectors of Ls, and number K of clusters
  %prompt = 'Method (1 = svd(Ls), 2 = svd(Bs), 3 = eig(Ls), 4 = svd(2*I - Ls); 0 quit:';
  %mm = input(prompt);
  mm=4;
  %prompt = 'Number of clusters (0 to end):';
  %K = input(prompt); 
  K=k;
  eigindex = 1;
  NNZ = 100;  % sets norm of Z to NNZ
  
  % Finds relaxed solution Z according to method mm
  Z1 = findZ(mm,m,K,U5,S5,V5,Ls,Bs,Dnhalf,0,NNZ,eigindex,0,0);
  
  %  To make the columns of Z Euclidean orthogonal as well as D orthogonal
  %  by finding R1 in O(K)
  %  [R1,SS,R0] = svd(Z,0);
  [R1,SSb] = eig(Z1'*Z1);
  Z2 = Z1*R1;
  
  %prompt = 'Initialize Z (1:basic, 2:rows sum to 1, 3:unit rows, 4:unit rows col ortho):';
  %initsw = input(prompt);
  initsw=1;
  %prompt = 'Method to find A, for min ||X - Z*A|| (1: R in O(k), 2: A in GL(k), 3: R*Lam):';
  %whichA = input(prompt); 
  whichA=findAmethod;
  
  % Tries to make the sum of the rows of Z1 and Z2 as close to rowsum as possible
  rowsum = 1;
  [Zr,~,fail] = rowsumto1(Z2,rowsum,0);
  if fail == 1
     fprintf('Failed in rowtosum1 for Z2 \n')
     Zr = Z2;
  end
  % Zr
  [Zr1,~,fail1] = rowsumto1(Z1,rowsum,0);
  if fail1 == 1
     fprintf('Failed in rowtosum1 for Z1 \n')
     Zr1 = Z1;
  end
  %  Normalizes the rows to have unit length
  Zu = normrows(Z2);
  Zu1 = normrows(Z1);
  % Tries to make the rows of unit length while preserving orthogonal columns
  [Zn,~,~,flag,~] = normcol(Z2);
  [Zn1,~,~,flag1,~] = normcol(Z1);

  if flag == 1
     fprintf('Failed in normcol for Z2 \n')
     Zn = Z2;
  end
  if flag1 == 1
     fprintf('Failed in normcol for Z1 \n')
     Zn1 = Z1;
  end
  if initsw == 2
     Zinit2 = Zr; Zinit1 = Zr1;
  else
      if initsw == 3
         Zinit2 = Zu; Zinit1 = Zu1;
      else
         if initsw == 4 
            Zinit2 = Zn; Zinit1 = Zn1;
         else    
            Zinit2 = Z2; Zinit1 = Z1; 
         end
      end
  end
  
  a0 = NNZ/sqrt(m); a = a0;     %  so that ||Xc||_F = NNZ
  if initsw == 3    
     a = 1;   a0= 1; % case where the rows of Zinit have unit length, so ||Zinit||_F = m 
  else
     NZ2 = sqrt(trace(Zinit2'*Zinit2));
     Zinit2 = (NNZ/NZ2)*Zinit2;   %  normalizes Zinit2 so that ||Zinit2||_F = NNZ 
     NZ1 = sqrt(trace(Zinit1'*Zinit1));
     Zinit1 = (NNZ/NZ1)*Zinit1;   %  normalizes Zinit1 so that ||Zinti1||_F = NNZ
  end
  
  step = 1;
  fprintf('step = %d \n',step)
  % showZ = 1;  %  to display the graph associated with Z
  % Finds the best initial rotation Rc and discrete solution Xc
  [Rc,sw,N1,N2,N3,N4,Xc] = initR2_v4(W,find2_X,Zinit2,Zinit1,R1,K,a0,a,0,0,0);
  if sw == 1  || sw == 2   %  the best initial solution is obtained with Zinit1
     Z = Zinit1; 
  else
     Z = Zinit2;
  end

  fprintf('*********************************************** \n')
  fprintf('Intitial R and X have been found; iterating to improve them \n')
  

  eX = sqrt(trace((Xc - Z*Rc)'*(Xc - Z*Rc)));

  if whichA == 2 
     [Rc,fail] = find_A(Xc,Z,show1,show2); 
     if fail == 1
        fprintf('A is singular in find_A  \n')
     end
  else
     if whichA == 3
        Rc = find_RLam(Xc,Z,show1,show2);
     else
        Rc = find_R(Xc,Z);   
     end
  end

  eRn = sqrt(trace((Xc - Z*Rc)'*(Xc - Z*Rc))); eR = eRn + 1;
  convlist = [];
  convlist(1) = eX; convlist(2) = eRn;  % convergence list
  Xc0 = Xc; Rc0 = Rc;
 
  %  finds best Xc and Rc
  while eRn < eR && step <= maxiter
     step = step + 1;
     fprintf('step = %d \n',step)
     [Xc, ~]  = find2_X(Z,Rc,a,0,0);  
     DiffXc = (Xc - Xc0)/a;

     NDXc = trace(DiffXc'*DiffXc)/2;  % Number of rows in Xc that changed
     if NDXc < tol
        nochangeinX = 1;
        fprintf('Xc did not change, NDXc = %d \n',NDXc)
     else
        nochangeinX = 0;
        fprintf('Xc changed, NDXc = %d \n',NDXc)
     end
     % Zn1 = Z*Rc; 
     eX = sqrt(trace((Xc - Z*Rc)'*(Xc - Z*Rc)));
     convlist(2*step-1) = eX;
     % eX1 = trace(Zn1'*Xc);
     if nochangeinX == 1
       % fprintf('Same Xc and Rc as in the previous step \n')   
        fprintf('Number of steps to reach a minimum Ncut = %d \n',step-1)
        convlist(2*step) = eRn;
        eRn = eR;  % stop the while loop
     else
        if whichA == 2
           [Rc,fail] = find_A(Xc,Z,0,0);
           if fail == 1
              fprintf('A is singular in find_A \n')              
           end
        else
           if whichA == 3
              Rc = find_RLam(Xc,Z,0,0);
           else
              Rc = find_R(Xc,Z);  
           end
        end
        eR = eRn;
        eRn = sqrt(trace((Xc - Z*Rc)'*(Xc - Z*Rc)));
        convlist(2*step) = eRn;
        if eR < eRn
           Xc = Xc0; Rc = Rc0;   % the last try was worse, use the previous one
           fprintf('Last try was worse than the previous step  \n')
           fprintf('Number of steps to reach a minimum Ncut = %d \n',step-1)
        else
           Xc0 = Xc; Rc0 = Rc;   % last try was better
        end
     end

     % eX2 = trace(Zn2'*Xc);
  end
  
  fprintf('*********************************************************** \n')
  fprintf('total number of steps = %d \n',step)
  XXc = (1/a)*Xc;           % to get a matrix of 0's and 1's
  [mval IDX] = max(XXc');
  IDX = IDX';
end

