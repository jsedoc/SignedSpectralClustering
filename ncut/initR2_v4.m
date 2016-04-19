%
% Function to compute the best initial rigid motion Rc in O(K)
% and initial discrete solution Xc
% This version tries in the following order:
% 1: Z1 and R3 = I;  sw = 1
% 2: Z1 and R2a from make_ortho_rows_v4; sw = 2
% 3: Z2 and R3 = I; sw = 3
% 4: Z2 and R2 from make_ortho_rows_v4; sw = 4
% Use improved make_ortho_row_v4
%  make_ortho_v5 picks the first row at random

function [Rc,sw,N1,N2,N3,N4,Xc] = initR2_v4(W,find2_X,Z2,Z1,R1,K,a0,a,show1,show2,showZ)
 % fprintf('a in initR2_v4 = %d \n',a)
 % fprintf('a0 in initR2_v4 = %d \n',a0)
 [R2,R2z] = make_ortho_rows_v4(Z2); 
 % [R2,R2z,~] = make_ortho_rows_v5(Z2);
  %  R2 has a zero column iff R2z = 1
  if show2 == 1
     R2
    % RR2 = R2'*R2;
    % RR2
  end

  [R2a,R2za] = make_ortho_rows_v4(Z1);
  
 % [R2a,R2za,~] = make_ortho_rows_v5(Z1);
  if show2 == 1
     R2a
    % RR2a = R2a'*R2a;
    % RR2a
  end
  R2b = R1*R2;
  R2c = R1'*R2a;
  if show2 == 1
     R2b
     R2c
  end
  DR = R2a - R2b;
  NDR = sqrt(trace(DR'*DR));
  if show2 == 1
     DR
     fprintf('NDR =||R1*R2 - R2a|| in initR2_v4 = %d \n',NDR)  
  end
  R3 = eye(K);
  % fprintf('a in initR2_v4 = %d \n',a)
  
  %  Finds closest Xc1 with Z1 and R3 = I
  if show2 == 1
     fprintf('Find Xc1 \n') 
  end
  [Xc1, RR1] = find2_X(Z1,R3,a0,show1,show2);
  % Xc1
  % Z1
  N1 = sqrt(trace((Xc1 - Z1*RR1)'*(Xc1 - Z1*RR1)));

  
  %  Finds closest Xc2 with Z1 and R2a
  if show2 == 1
     fprintf('Find Xc2 \n') 
  end
  [Xc2, RR2] = find2_X(Z1,R2a,a0,show1,show2);
  % Xc2
  % Z1
  N2 = sqrt(trace((Xc2 - Z1*R2a*RR2)'*(Xc2 - Z1*R2a*RR2)));

  
  %  Finds closest Xc3 with Z2 = Z1*R1 and R3 = I
  if show2 == 1
     fprintf('Find Xc3 \n') 
  end
  [Xc3, RR3] = find2_X(Z2,R3,a,show1,show2);
  N3 = sqrt(trace((Xc3 - Z2*RR3)'*(Xc3 - Z2*RR3)));

  
  %  Finds closest Xc4 with Z2 and R2
  if show2 == 1
     fprintf('Find Xc4 \n') 
  end
  [Xc4, RR4] = find2_X(Z2,R2,a,show1,show2);
  % RR4 is applied second because it flips the signs of columns of Z2*R2
  N4 = sqrt(trace((Xc4 - Z2*R2*RR4)'*(Xc4 - Z2*R2*RR4)));

  % NOTE: type 2 - 4 ZZ^T vs pick k most ortho 
  DiffZR = Z2*R2*RR4 - Z1*R2a*RR2;
  NDiff24 = sqrt(trace(DiffZR'*DiffZR));
  if show2 == 1
      fprintf('NDiff24 in initR2_v4 = %d \n',NDiff24)
      fprintf('N1 = %d \n',N1)
      fprintf('N2 = %d \n',N2)
      fprintf('N3 = %d \n',N3)
      fprintf('N4 = %d \n',N4)
  end
  [minN,sw] = min([N1 N2 N3 N4]);
  if show2 == 1
      fprintf('MinN in initR2_v4 = %d \n',minN)
      fprintf('sw = %d \n',sw)
  end
  if sw == 1
      Rc = RR1;  Xc = Xc1; 
  else
      if sw == 2
         Rc = R2a*RR2; Xc = Xc2;          
      else
         if sw == 3
            Rc =  RR3; Xc = Xc3;    
         else
            Rc = R2*RR4; Xc = Xc4;
         end
      end
  end  
  if show1 == 1
     fprintf('Intial R and X in initR2_v4 \n')
     Rc
     Xc
  end
  if R2z == 1 || R2za == 1
     if show1 == 1  
       fprintf('R2 has some zero column \n')  
     end
     sw = 1;
  end  
end