%
% Normalizes the rows of Z2
%

function Zu = normrows(Z2, fast)
% Check number of inputs.
if nargin > 2
    error('myfuns:normrow:TooManyInputs', ...
        'requires at most 2 inputs');
end

if nargin == 1
    fast=1;
end

% NOTE: element wise - this is N^2
if fast
    Zu = bsxfun(@rdivide, Z2, sqrt(sum(Z2.^2, 2)));
else
 m = size(Z2,1);
 Z2rowN = Z2*Z2'; Z2rowinv = zeros(1,m);
 for i = 1:m 
     Z2rowinv(i) = Z2rowN(i,i)^(-1/2);
 end
 Z2rowi = diag(Z2rowinv);
 Zu = Z2rowi*Z2;
  % NZu = Zu*Zu';
end
end

