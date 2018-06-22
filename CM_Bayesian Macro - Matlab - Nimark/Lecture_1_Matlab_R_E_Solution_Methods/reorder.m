function [s,t,q,z] = reorder(s,t,q,z)

%
% Takes U.T. matrices S, T, orthonormal matrices Q,Z, rearranges them
% so that abs(T(i,i)/S(i,i)) are in ascending order
% while preserving U.T. and orthonormal properties and Q'AZ' and
% Q'BZ'.
%


n = size(s,1);
i = 1;
while i<=n-1;
   if abs(t(i,i)*s(i+1,i+1))>abs(s(i,i)*t(i+1,i+1));    
      [s,t,q,z] = qzswitch(i,s,t,q,z);
      if ~(i==1);i = i-2;end
   end
   i=i+1;
end

   



