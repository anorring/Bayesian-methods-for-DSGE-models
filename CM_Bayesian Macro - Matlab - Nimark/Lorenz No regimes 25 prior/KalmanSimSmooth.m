function [Xss]=KalmanSimSmooth(theta,A,C,Z,Atilda,Btilda,SigJ,dimX,T,yc,dimx,NonZeroBond,maxsurv,H,dimu,muP)
global FullInfoAllSurveys

%  Kalman simulation smoother
%  Adapted from Durbin and |Koopman (2002) by K Nimark

%     Xt = A*Xt-1 +C*ut
%
%     Zt=  D1*Xt + R*ut


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First demean the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimobs=length(yc);
CC=(C*C');

P0=CC;
P0st=P0;
tolerance=1d-6;
iter=0;
diff=1;

while (diff>tolerance)   
    iter=iter+1;
    P0 = A*P0st*A' + CC;
    diff=max(max(abs(P0st-P0)./(1 + abs(P0))));
    P0st=P0;    
    
    if (iter>25000)
        diff=0;
        disp('more than 25000 iterations are needed to compute unconditional variance')
    end    
    
end

Lambda0=zeros(dimX,1);
Lambda0(1:dimx,1)=theta(15:17);

muX=eye(dimX)/(eye(dimX)-A)*muP;

ZZ(1:dimobs,:)=Z(1:dimobs,:)-kron(Atilda(yc),ones(1,T)) - kron(Btilda(yc,:)*muX,ones(1,T));
ZZ(dimobs+1:dimobs+maxsurv,:)=Z(dimobs+1:dimobs+maxsurv,:)-Atilda(1)-Btilda(1,:)*muX;
Dsurv=Btilda(1,:)*A;
survsig=theta(end);

% ZZ(1:dimobs,:)=Z(1:dimobs,:)-kron(Atilda(yc),ones(1,T)) - kron(Btilda(yc,:)/(eye(dimX)-A)*CC*Lambda0,ones(1,T));
% ZZ(dimobs+1:dimobs+maxsurv,:)=Z(dimobs+1:dimobs+maxsurv,:)-Atilda(1)-Btilda(1,:)/(eye(dimX)-A)*CC*Lambda0;
% Dsurv=Btilda(1,:)*A;
% survsig=theta(end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define a ancilliary variables, predefine matrices etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimS=dimu+maxsurv;   %dimension of shocks;
dimZ=dimobs+maxsurv;   %dimension of observables;

Xhat=zeros(dimX,T+1);
PP0=zeros(dimX,dimX,T);
PP1=zeros(dimX,dimX,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%draw from (unconditional) state distibution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xplus=zeros(dimX,T+1);
%Zplus=NaN(dimZ,T+1);
Xplus(:,1)=chol(P0+(1e-14)*eye(dimX))'*randn(dimX,1);
shocks=randn(dimS,T);
for t=1:T;
    
    %RRR = [zeros(dimobs,dimu) zeros(dimobs,NonZeroBond(t));
    %zeros(NonZeroBond(t),dimu) eye(NonZeroBond(t))*(Dsurv*SigJ*Dsurv')^0.5 ];

    %RRR(2:dimobs,dimx+yc(2:end)-1)=theta(13)*eye(dimobs-1);   %Since the short rate is observed we look at yc(2:end)

    NRR=[C zeros(dimX,NonZeroBond(t))];
    
    Xplus(:,t+1) = A*Xplus(:,t) + NRR*shocks(1:dimu+NonZeroBond(t),t);
    %Zplus(:,t+1) = RRR*shocks(1:dimu+NonZeroBond(t),t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct new observables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zstar=ZZ;%-Zplus(:,2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forward recursion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt=1:T
    
    D1=[Btilda(yc,:);
       ones(NonZeroBond(tt),1)*Dsurv*H;];

    if (FullInfoAllSurveys)
            R = [zeros(dimobs,dimu) zeros(dimobs,NonZeroBond(tt));
               zeros(NonZeroBond(tt),dimu) eye(NonZeroBond(tt))*survsig];
    else
        R = [zeros(dimobs,dimu) zeros(dimobs,NonZeroBond(tt));
               zeros(NonZeroBond(tt),dimu) eye(NonZeroBond(tt))*(Dsurv*SigJ*Dsurv')^0.5 ];
    end
      
    %R(2:dimobs,dimx+yc(2:end)-1)=theta(13)*eye(dimobs-1);   %Since the short rate is observed we look at yc(2:end)
    R(1:dimobs,(dimx+1):(dimx+dimobs))=theta(13)*eye(dimobs);   %Since the short rate is observed we look at yc(2:end)
    
    Ztilde=Zstar(1:dimobs+NonZeroBond(tt),tt)-D1*A*Xhat(:,tt);
 
    NRR=[C zeros(dimX,NonZeroBond(tt))];
    Omega=(D1*A)*P0*(D1*A)'+(D1*NRR+R)*(D1*NRR+R)';
    Omegainv=eye(dimobs+NonZeroBond(tt))/Omega;
    K=(A*P0*(D1*A)'+CC*D1'+NRR*R')*Omegainv;
    Xhat(:,tt+1)=A*Xhat(:,tt)+K*Ztilde;
    P1=A*P0*A'+CC;
    P0=P1-K*Omega*K';
    P1=A*P0*A'+CC;
    PP0(:,:,tt)=P0;
    PP1(:,:,tt)=P1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct what is needed for last step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xtilde=A*Xhat(:,1:end-1);
Xhat=Xhat(:,2:end);

Xsm=Xhat*0;
Xsm(:,T)=Xhat(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%backward recursion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tt=T:-1:2
    J=PP0(:,:,tt-1)*A'*inv(PP1(:,:,tt-1)+eye(dimX)*1e-6);
    Xsm(:,tt-1)=Xhat(:,tt-1)+J*(Xsm(:,tt)-Xtilde(:,tt));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add up and spit out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xss=Xsm+Xplus(:,2:end)+kron(muX,ones(1,T));



% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %draw from (unconditional) state distibution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xplus=zeros(dimX,T+1);
% Zplus=zeros(dimZ,T);
% Xplus(:,1)=chol(P0+(1e-2)*eye(dimX))'*randn(dimX,1);
% shocks=randn(dimS,T);
% for tt=2:T+1;
%     Xplus(:,tt)=A*Xplus(:,tt-1) + C*shocks(:,tt-1);
%     Zplus(:,tt) = R*shocks(:,tt-1);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Construct new observables
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zstar=Z;%-Zplus(:,2:end);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % forward recursion
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for tt=1:T
%     
%     
%     Ztilde=Zstar(:,tt)-D1*A*Xhat(:,tt);
%     Omega=(D1*A)*P0*(D1*A)'+(D1*C+R)*(D1*C+R)';
%     Omegainv=eye(dimZ)/Omega;
%     K=(A*P0*(D1*A)'+C*C'*D1'+C*R')*Omegainv;
%     Xhat(:,tt+1)=A*Xhat(:,tt)+K*Ztilde;
%     P1=A*P0*A'+C*C';
%     P0=P1-K*Omega*K';
%     P1=A*P0*A'+C*C';
%     PP0(:,:,tt)=P0;
%     PP1(:,:,tt)=P1;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %construct what is needed for last step
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Xtilde=A*Xhat(:,1:end-1);
% Xhat=Xhat(:,2:end);
% 
% Xsm=Xhat(:,2:end)*0;
% Xsm(:,T)=Xhat(:,end);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %backward recursion
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for tt=T:-1:2
%     J=PP0(:,:,tt-1)*A'*inv(PP1(:,:,tt-1)+eye(dimX)*1e-14);
%     Xsm(:,tt-1)=Xhat(:,tt-1)+J*(Xsm(:,tt)-Xtilde(:,tt));
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %add up and spit out
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xss=Xsm+Xplus(:,2:end);
