function [K,P,p,Xtt]=kalman_filter(A,C,D,R,Z)
%Computes Kalman filter for state space system
% 	X[t] = AX[t-1] + Cu[t]
% 
%   Z[t] = DX[t] + Ru[t]
%
% Outputs are the Kalman gain K in 
% X[t|t] = X[t|t-1]+ K(Z[t] - DX[t|t-1] )
%
% P is the steady state prior error covariance matrix
%
%  P = E(X[t]-X[t|t-1])(X[t]-X[t|t-1])'
% 
%  p is the steady state posterior error covariance matrix
%
%  p = E(X[t]-X[t|t])(X[t]-X[t|t])'
% and filtereed estiamte of state X[t|t]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimX=length(A(:,1));
dimZ=length(D(:,1));
periods=length(Z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Kalman filter equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=zeros(dimX,dimX,periods+1);
K=zeros(dimX,dimZ,periods);
Xtt=zeros(dimX,1);
P(:,:,1)=C*C';
P(:,:,1)=C*C'/(1-A^2);
for t=1:periods;
    P(:,:,t+1)=A*(P(:,:,t)-(P(:,:,t)*D'/(D*P(:,:,t)*D'+R*R'))*D*P(:,:,t))*A'+C*C';    
    K(:,:,t)= P(:,:,t+1)*D'/(D*P(:,:,t+1)*D'+R*R');
    Xtt(:,t+1)=A*Xtt(:,t)+K(:,:,t)*(Z(:,t)-D*A*Xtt(:,t));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute steady state Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1st=C*C';
diff=1;
iter=1;
tol=1e-6;
maxiter=1e4;
while diff>= tol && iter <= maxiter
    P1=A*(P1st-(P1st*D'/(D*P1st*D'+R*R'))*D*P1st)*A'+C*C';    
    diff=max(max(abs(P1-P1st)));
    iter=iter+1;
    P1st=P1;
end
K=P1*D'/(D*P1*D'+R*R');
P=P1;
p=(P1st-(P1st*D'/(D*P1st*D'+R*R'))*D*P1st);
 