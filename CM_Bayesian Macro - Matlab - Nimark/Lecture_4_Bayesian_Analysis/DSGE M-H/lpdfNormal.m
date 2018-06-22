function  LogDens = lpdfNormal(x,a,b)

%--------------------------------------------------------------------------
% Compute the log of the N(mu,sigma^2) density
%
%  x: Evaluated at  
%  a: Mean
%  b: Standard deviation
%      
%--------------------------------------------------------------------------

LogDens = -log(b) -0.5*log(2*pi) -0.5*((x-a)./b).^2;
 
