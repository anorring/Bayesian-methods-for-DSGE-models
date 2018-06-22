										
function[K,S]=doubleo(A,C,Q,R)
%DOUBLEo.M
%function[K,S]=doubleo(A,C,Q,R)
%  This program uses the "doubling algorithm" to solve the
%  Riccati matrix difference equations associated with the
%  Kalman filter.  A is nxn, C is kxn, Q is nxn, R is kxk.
%  The program returns the gain K and the stationary covariance
%  matrix of the one-step ahead errors in forecasting the state.
%
%  The program creates the Kalman filter for the following system:
%
%       x(t+1) = A * x(t) + e(t+1)
%
%         y(t) = C * x(t) + v(t)
%
%  where E e(t+1)*e(t+1)' =  Q, and E v(t)*v(t)' = R, and v(s) is orthogonal
%  to e(t) for all t and s.  The program creates the observer system
%
%        xx(t+1) = A * xx(t) + K * a(t)
%           y(t) = C * xx(t) + a(t),
%
%  where K is the Kalman gain ,S = E (x(t) - xx(t))*(x(t) - xx(t))', and
%  a(t) = y(t) - E[y(t)| y(t-1), y(t-2), ... ], and xx(t)=E[x(t)|y(t-1),...].
%  NOTE:  By using DUALITY, control problems can also be solved.
a0=A';
b0=C'*(R\C);
g0=Q;
tol=1e-15;
dd=1;
ss=max(size(A));
v=eye(ss);
while dd>tol
a1=a0 *((v+b0*g0)\a0);
b1=b0+a0*((v+b0*g0)\(b0*a0'));
g1=g0+a0'*g0*((v+b0*g0)\a0);
k1=A*g1*C'/(C*g1*C'+R);
k0=A*g0*C'/(C*g0*C'+R);
dd=max(max(abs(k1-k0)));
a0=a1;
b0=b1;
g0=g1;
end
K=k1;S=g1;
