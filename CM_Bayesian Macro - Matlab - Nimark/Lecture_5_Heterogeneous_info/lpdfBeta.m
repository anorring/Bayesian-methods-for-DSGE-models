function LogDens = lpdfBeta(x,a,b);

%--------------------------------------------------------------------------
% Compute the logBeta(a,b) density fcn
%
%      x:     The density is evaluated at x  
%      a:     First (converted) parameter
%      b:     Second (converted) parameter
%      
%--------------------------------------------------------------------------


 if x>=0 & x<=1
     LogDens = -betaln(a,b) + (a-1).*log(x) + (b-1).*log(1-x);
 else
     LogDens=-inf;
 end