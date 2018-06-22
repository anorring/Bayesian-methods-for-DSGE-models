function LogDens = lpdfInvGamma(x,a,b)

%--------------------------------------------------------------------------
% The log height of the InvertedGamma(a,b) density
%
%      x:  The density is evaluated at  
%      a:  First parameter
%      b:  Second parameter
%       
%--------------------------------------------------------------------------

LogDens = log(2)-gammaln(b/2)+(b/2).*log(b*a^2/2) -((b+1)/2).*log(x.^2)-b*a^2./(2*x.^2);

