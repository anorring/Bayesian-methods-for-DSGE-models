function [X]=sim_smooth(A,C,D,R,Z)

%  Kalman simulation smoother
%  Adapted from Durbin and |Koopman (2002) by K Nimark

%     Xt = A*Xt-1 +C*ut
%
%     Zt=  D*Xt + R*ut


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define a ancilliary variables, predefine matrices etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=length(Z);
dimX=length(A);        %dimension of state
dimS=length(R(1,:));   %dimension of shocks;
dimZ=length(D(:,1));   %dimension of observables;

CC=C*C';

if abs(max(eig(A)))>=1;
P0=CC*1000;
else
P0=reshape(inv(eye(dimX*dimX)-kron(A,A))*CC(:),dimX,dimX);
end
Xhat=zeros(dimX,T+1);
PP0=zeros(dimX,dimX,T);
PP1=zeros(dimX,dimX,T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%draw from (unconditional) state distibution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xplus=zeros(dimX,T+1);
Zplus=zeros(dimZ,T);
Xplus(:,1)=chol(P0+1e-6*eye(dimX))'*randn(dimX,1);
shocks=randn(dimS,T);
for tt=2:T+1;
    Xplus(:,tt)=A*Xplus(:,tt-1) + C*shocks(:,tt-1);
    Zplus(:,tt) = D*Xplus(:,tt) + R*shocks(:,tt-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct new observables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zstar=Z;%-Zplus(:,2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forward recursion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt=1:T
    Ztilde=Zstar(:,tt)-D*A*Xhat(:,tt);
    Omega=(D*A)*P0*(D*A)'+(D*C+R)*(D*C+R)';
    Omegainv=eye(dimZ)/Omega;
    K=(A*P0*(D*A)'+C*C'*D'+C*R')*Omegainv;
    Xhat(:,tt+1)=A*Xhat(:,tt)+K*Ztilde;
    P1=A*P0*A'+C*C';
    P0=P1-K*Omega*K';
    P1=A*P0*A'+C*C';
    PP0(:,:,tt)=P0;
    PP1(:,:,tt)=P1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct what is needed for last step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xtilde=A*Xhat;
Xsm=Xhat*0;
Xsm(:,T)=Xhat(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%backward recursion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tt=T:-1:1
    J=PP0(:,:,tt)*A'*inv(PP1(:,:,tt)+eye(dimX)*1e-10);
    Xsm(:,tt)=Xhat(:,tt)+J*(Xsm(:,tt+1)-Xtilde(:,tt+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add up and spit out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=Xsm+Xplus;
