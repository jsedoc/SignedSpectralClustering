%
% To compute relaxed solution Z according to method mm
%

function [Z,quit] = findZ(mm,m,K,U5,S5,V5,Ls,Bs,Dnhalf,quit,NNZ,eigindex,show1,show2)
    Z = 0;
if mm == 1
     if (issparse(Ls))
        [Un,Sn,Vn] = svds(Ls,size(Ls,1));        % SVD of Normalized Laplacian
     else
        [Un,Sn,Vn] = svd(Ls);        % SVD of Normalized Laplacian
     end

    if show2 == 1
       Vn
       Un
       Sn
    end
    U = Un;       %  relaxed solution using an SVD of the Laplacian Ls
    Ls2b = Un*Sn*Un';
    DifL2s = sqrt(trace((Ls2b - Ls)'*(Ls2b - Ls)));
    % Un
  else
    if mm == 2  
%      [U1, S1] = svd(Bs,'econ');   % SVD of generalized incidence matrix
      [U1, S1, V1] = svd(Bs);       % SVD of generalized incidence matrix
      if show2 == 1
          V1
          U1
          S1
      end
       U = U1;       %  relaxed solution using an SVD of the incidence matrix  Bs
       S2 = S1*S1';
       % S2 = S1(:,1:min(m,n))^2;
       Ls3 = U1*S2*U1';  % normalized  Laplacian from SVD of normalized incidence matrix
       DifLs3 = sqrt(trace((Ls3 - Ls)'*(Ls3 - Ls)));
       % U1
    else
        if mm == 3
           [v, e] = eig(Ls);  %   eigenvectors of Ls
           if show2 == 1
               v
               e
           end
      %    U = v(:,1:K);
           U = v(:,eigindex:eigindex+K-1);
           Ls1b = v*e*v';
           DifLs1 = sqrt(trace((Ls1b - Ls)'*(Ls1b - Ls)));
           % v
        else
            if mm == 4
              %  ILs = 4*eye(m) - Ls;
              %  ILs = inv(eye(m) + Ls);
              %  [U4,S4] = eigs(ILs,K,'la');
              %  [U5,S5,V5] = svd(ILs);
              %  [U5,S5] = eig(ILs);
                if show2 ==1
                    V5
                    U5
                    S5
                end
              %  U = U5(:,1:K);
                 U = U5(:,eigindex:eigindex+K-1);
              %  Ls4 = U5*S5*U5';
                Ls4 = U5*S5*U5';
                DifLs4 = sqrt(trace((Ls4 + Ls -4*eye(m))'*(Ls4 + Ls - 4*eye(m))));
                % U
            else
                quit = 1; return
            end
        end
    end
  end
  if mm == 3 || mm == 4
      Y = U;
  else
      % Y = U(:,m-K+1:m);
        Y = U(:,m-K+2-eigindex:m+1-eigindex);
  end
  % Y
  Z = Dnhalf*Y;     % relaxed solution
  NZ = sqrt(trace(Z'*Z)); 
  Z = (NNZ/NZ)*Z;   %  normalizes Z so that ||Z||_F = NNZ 
  if show1 == 1
     fprintf('Solution of the relaxed problem Z \n')
     Z
  end
  if show2 == 1   % shows error between Ls and the other Laplacians
     if mm == 1
        fprintf('||Un*Sn*Un^T - Ls|| = %d \n',DifL2s)
     else
        if mm == 2
           fprintf('||U*S^2*U^T - Ls|| = %d \n',DifLs3)
        else
           if mm == 3
              fprintf('||v*e*v^T - Ls|| = %d \n',DifLs1)
           else
              fprintf('||U5*S5*U5^T - Ls|| = %d \n',DifLs4)
           end
        end
     end 
  end
  quit = 0;
end

