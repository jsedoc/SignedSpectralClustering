%
% Function implementing the normalized cut method for K clusters
% This version picks the best initialization for R
% and accomodates negative weights by using absolute values
% in the degree matrix. In this version, displaying the various
% graphs is a separate function
% rw = 1 means random perturbation added to graph nodes
% Given the intial relaxed solution Z1, R1 in O(k) is found
% from the svd  R1*Sigma*R1^T of Z1^T*Z1. Then, Z2 = Z1*R1
% is both Euclidean otho and D-ortho, where D is the (signed) degree
% matrix.
%  This version uses three different initialization methods for R
%  1:  uses Z2 = Z1*R1 and intial R = I
%  2:  uses Z2 = Z1*R1 and R2 obtained from Z2 using make_ortho_v4
%  3:  uses Z1 and R2a obtained from Z1 using make_ortho_v4. In principle
%      R2a = R1*R2, but due to numerical errors, this may not be the case 
%  4:  uses Z1 and R = I
%  Also, tries different versions of Z2
%  initsw = 1: basic method using Z2
%  initsw = 2: tries to make rows sum to 1; uses Zr 
%  initsw = 3: normalizes the rows of Z2; uses Zn
%  initsw = 4: tries to make the rows of unit length by recaling the
%              columns; calls normcal. It is fails, uses Z2.
%  specifies index of first eigenvector as eigindex
%  Best method to compute eigenvalues of Ls: method 4, use 
%  4*I - Ls.
%  This version uses two other methods to find R
%  2: find_A looks for an invertible
%  matrix such that ||X - Z*A||_F is minimized
%  3: find_RLam looks for a matrix in O(K) and a diagonal matrix Lam
%  such that ||X - Z*R*Lam||_F is minimized
%  showZ is used in showclust_v4 to display the graph associated with Z
%  displays iff showZ = 1

function [L,Ls,Bs,Z1,Z2,Zinit1,Zinit2,R1,Xc,Rc] = ncutK_v5(W,nx,ny,show1,show2,showZ)
% sw = swicth to decide how to initialize Rc
m = size(W,1);   tol = 10^(-14); % tolerance to decide when an eigenvalue is 0
minval = min(min(W));   % finds smallest entry in W
if minval < 0
     fprintf('The graph is a signed graph  \n')
else
     fprintf('The graph is an unsigned graph  \n')
end
Dv = sum(abs(W));       % vector of degrees
D = diag(Dv);           %  computes degree matrix of absolute values
figure(1)
%colormap(gray)
imagesc(50*W)
vertlist = drawgraph_v2(nx,ny);
show_the_graph(W,nx,ny,vertlist)
degv = D*ones(m,1); d = ones(m,1)'*degv; % vector of degrees and total volume
%W
% Dhalfv = Dv.^(1/2);  Dhalf = diag(Dhalfv);
Dnhalfv = Dv.^(-1/2); Dnhalf = diag(Dnhalfv);
% Dnhalf
L = D - W;             % standard signed Laplacian
% L
Ls = Dnhalf*L*Dnhalf;  % normalized  Laplacian
B = incidmat_v2(W);    % generalized incidence matrix of W
% B
n = size(B,2);
Bs = Dnhalf*B;         % normalized incidence matrix
% DifLs = Bs*Bs' - Ls;
% DifLs

quit = 0;
while quit ~= 1
  fprintf('Number of nodes = %d \n',m)
  fprintf('Number of edges = %d \n',n)
  [U5,S5,V5] = svd(2*eye(m) - Ls);        % SVD of 2*I - Ls
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
  prompt = 'Method (1 = svd(Ls), 2 = svd(Bs), 3 = eig(Ls), 4 = svd(2*I - Ls); 0 quit:';
  mm = input(prompt);
  prompt = 'Number of clusters (0 to end):';
  K = input(prompt); 
  if K < 2 || K > m 
     quit = 1; break
  end
  prompt = 'first eigenvector index (standard = 1):';
  eigindex = input(prompt);
  if eigindex < 1  || eigindex > m - K + 1
     eigindex = 1;
  end
  NNZ = 100;  % sets norm of Z to NNZ
  
  % Finds relaxed solution Z according to method mm
  [Z1,quit] = findZ(mm,m,K,U5,S5,V5,Ls,Bs,Dnhalf,quit,NNZ,eigindex,show1,show2);
  if show2 == 1
     fprintf('relaxed solution Z1  \n')
     Z1
  end
  if quit == 1 
     break
  end    
  %  To make the columns of Z Euclidean orthogonal as well as D orthogonal
  %  by finding R1 in O(K)
  %  [R1,SS,R0] = svd(Z,0);
  [R1,SSb] = eig(Z1'*Z1);
  if show2 == 1
     fprintf('rigid motion R1 to make Z2 = Z1*R1 ortho and D-ortho  \n')
     R1 
  end
  Z2 = Z1*R1;
  if show2 == 1
     fprintf('Z2 = Z1*R1  \n')
     Z2
  end
  prompt = 'Initialize Z (1:basic, 2:rows sum to 1, 3:unit rows, 4:unit rows col ortho):';
  initsw = input(prompt);
  prompt = 'Method to find A, for min ||X - Z*A|| (1: R in O(k), 2: A in GL(k), 3: R*Lam):';
  whichA = input(prompt); 
  % Tries to make the sum of the rows of Z1 and Z2 as close to rowsum as possible
  rowsum = 1;
  [Zr,~,fail] = rowsumto1(Z2,rowsum,show2);
  if fail == 1
     fprintf('Failed in rowtosum1 for Z2 \n')
     Zr = Z2;
  end
  % Zr
  [Zr1,~,fail1] = rowsumto1(Z1,rowsum,show2);
  if fail1 == 1
     fprintf('Failed in rowtosum1 for Z1 \n')
     Zr1 = Z1;
  end
  if show1 == 2
     SZ = Zr*ones(K,1);
     SZ1 = Zr1*ones(K,1);
     SZ
     SZ1
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
  [Rc,sw,N1,N2,N3,N4,Xc] = initR2_v4(W,Zinit2,Zinit1,R1,K,a0,a,show1,show2,showZ);
  if sw == 1  || sw == 2   %  the best initial solution is obtained with Zinit1
     Z = Zinit1; 
  else
     Z = Zinit2;
  end
  if show2 == 1
     fprintf('Zinit1, Zinit2, Z afer initialization step, initsw = %d \n',initsw)
     Zinit1
     Zinit2
     Z
  end
  fprintf('*********************************************** \n')
  fprintf('Intitial R and X have been found; iterating to improve them \n')
  
  maxiter = 50;    
  % [Xcb, ~]  = find2_X(Z,Rc,a,show1,show2); 
  % DiffXc = Xc - Xcb;
  % DiffXc
  eX = sqrt(trace((Xc - Z*Rc)'*(Xc - Z*Rc)));
  if K > 2
     fprintf('Initial Xc and Rc  \n')
     figure
     text = 'Initial Xc and Rc';
     showclust_v4(W,Z,Xc,Rc,text,showZ);    % initial Xc and Rc
  end
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
  if K > 2
     fprintf('Xc and Rc after step %d  \n',step)
     figure
     text = 'Xc and Rc after step 1';
     showclust_v4(W,Z,Xc,Rc,text,showZ);   % Xc and Rc after step 1
  end
  eRn = sqrt(trace((Xc - Z*Rc)'*(Xc - Z*Rc))); eR = eRn + 1;
  convlist = [];
  convlist(1) = eX; convlist(2) = eRn;  % convergence list
  Xc0 = Xc; Rc0 = Rc;
 
  %  finds best Xc and Rc
  while eRn < eR && step <= maxiter
     step = step + 1;
     fprintf('step = %d \n',step)
     [Xc, ~]  = find2_X(Z,Rc,a,show1,show2);  
     DiffXc = (Xc - Xc0)/a;
     if show1 == 1
        DiffXc
     end
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
           [Rc,fail] = find_A(Xc,Z,show1,show2);
           if fail == 1
              fprintf('A is singular in find_A \n')              
           end
        else
           if whichA == 3
              Rc = find_RLam(Xc,Z,show1,show2);
           else
              Rc = find_R(Xc,Z,show1);  
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
           if K > 2
              figure
              text = {'Xc and Rc after step '}; strstep = {int_to_char(step)};
              text = strcat(text,strstep);
              showclust_v4(W,Z,Xc,Rc,text,showZ);
           end
        end
     end
     if show1 == 1
        Xc
        Zn2 = Z*Rc;
        Zn2
     end
     % eX2 = trace(Zn2'*Xc);
  end
  ZR = Z*Rc;
  
  fprintf('*********************************************************** \n')
  fprintf('total number of steps = %d \n',step)
  convlist
  numblocks = zeros(1,K);
  XXc = (1/a)*Xc;           % to get a matrix of 0's and 1's
  if show1 == 1
     fprintf('Solution of the discrete problem XXc \n')
     XXc
  end
  %  computes the number of blocks
  Nb = round(XXc'*XXc);
  for k = 1:K
      numblocks(1,k) = Nb(k,k);
  end
  % computes the value of the normalized cut
  Ncut = 0;
  for k = 1:K
      Ncut = Ncut + (XXc(:,k)'*Ls*XXc(:,k))/(XXc(:,k)'*D*XXc(:,k));
  end
  % fprintf('Number of nodes = %d \n',m)
  if show2 == 1
      fprintf('N1 (Z1 and R = I) = %d \n',N1)
      fprintf('N2 (Z1 and R = R2a) = %d \n',N2)
      fprintf('N3 (Z2 and R = I) = %d \n',N3)
      fprintf('N4 (Z2 and R = R2) = %d \n',N4)
  end
  fprintf('Initialization method for rotation R = %d \n',sw)
  fprintf('Intialization method for Z = %d \n',initsw)
  fprintf('Number of clusters = %d \n',K)
  numblocks
  fprintf('Ncut = %d \n',Ncut)
  prompt = 'Method to display partition  (1 with edges, 2 only nodes):';
  sw2 = input(prompt); 
  prompt = 'Add random perturbation to nodes to display graphs (yes = 1):';
  rw = input(prompt);
  show_graphs(W,K,XXc,numblocks,nx,ny,sw2,rw)
end
end

