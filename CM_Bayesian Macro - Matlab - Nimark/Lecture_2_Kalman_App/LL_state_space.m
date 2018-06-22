function [LL]=LL_state_space(A,C,D,R,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes log likelihood for state space system
% 	X[t] = AX[t-1] + Cu[t]
%
%   Z[t] = DX[t] + Ru[t]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial values for Kalman filter etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=length(Z);
LL=0;
if max(abs(eig(A)))<1
    P0=dlyap(A,C*C');
else
    P0=0;
    for j=1:100;
        P0=A*P0*A'+C*C';
    end
end
dimZ=length(D(:,1));
Xfilt=A(:,1)*0; %initial value for the filtered state.
%Compute recursive likleihood using the Kalman filter
for tt=1:T
    Ztilde=Z(:,tt)-D*A*Xfilt;
    Omega=(D*A)*P0*(D*A)'+(D*C+R)*(D*C+R)';
    Omegainv=eye(dimZ)/Omega;
    K=(A*P0*(D*A)'+C*C'*D'+C*R')*Omegainv;
    Xfilt=A*Xfilt+K*Ztilde;
    P1=A*P0*A'+C*C';
    P0=P1-K*Omega*K';
   
    LL = LL - 0.5*dimZ*log(2*pi) - 0.5*(log(det(Omega)) + Ztilde'*Omegainv*Ztilde);
end
