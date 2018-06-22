%Computes the Kalman filter for systems with lagged observables
function [K,P1,P0]=Kalman(A,C,D1,R)
% R=[R1 R2];
tol=1e-8;
maxiter=1000;
diff=1;iter=1;
P0=C*C';
P1=A*P0*A'+C*C';
while diff>= tol && iter <= maxiter 
   L=(D1*A)*P0*(D1*A)'+(D1*C+R)*(D1*C+R)';
   K=(A*P0*(D1*A)'+C*C'*D1'+C*R')/(L);
   P0=P1-K*L*K';
   P1st=A*P0*A'+C*C';
   diff=max(max(abs(P1-P1st)));
   iter=iter+1;
   P1=P1st;
end

P=P1;

    