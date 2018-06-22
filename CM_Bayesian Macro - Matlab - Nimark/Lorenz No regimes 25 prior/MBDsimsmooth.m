function [X]=MBDsimsmooth(M,N,a,b,e1,bindim,dimX,jlead1,jlead0,SigJ,theta,Z,ST,T,CPImean,NGDPmean,NonZeroCPI,NonZeroNGDP,maxsurv,H)
%  Kalman simulation smoother
%  Adapted from Durbin and |Koopman (2002) by K Nimark

%     Xt = A*Xt-1 +C*ut
%
%     Zt=  D1*Xt + D2*Xt-1 + R*ut


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define a ancilliary variables, predefine matrices etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xhat=zeros(dimX,T+1);
Xtilde=zeros(dimX,T+1);
PP0=zeros(dimX,dimX,T);
PP1=zeros(dimX,dimX,T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%draw from (unconditional) state distibution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pp=M(:,:,2)*0;
for j=1:10;
    pp=M(:,:,2)*pp*M(:,:,2)'+N(:,:,2)*N(:,:,2)';
end
Xplus=zeros(dimX,T+1);
Xplus(:,1)=chol(pp+1e-12*eye(dimX))'*randn(dimX,1);
dimS=length(N(1,:,1));
shocks=randn(dimS,T);

for tt=2:T;
    s=ST(tt-1:tt-1+bindim-2)';
    j=binvec2dec(s)+1;
    Xplus(:,tt)=M(:,:,j)*Xplus(:,tt-1) + N(:,:,j)*shocks(:,tt-1);
    %     Zplus(:,tt) = D1*Xplus(:,tt)+ D2*Xplus(:,tt-1) + R*shocks(:,tt-1);
end
   stick=theta(15); %Calvo parameter
beta=theta(16); %discount rate
delta=theta(10); %labour supply curvature

lambda=(1-stick)*(1-stick*beta)/stick;
Rlag=theta(12)*Z(2,1);
    Gr=[-1;-lambda+lambda*delta;]*theta(12)/(1-theta(12));
    Gu=[-1;-lambda+lambda*delta;]/(1-theta(12));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forward recursion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:T
    s=ST(t:t+bindim-1)';
    j=binvec2dec(s)+1;
    
    Fpi=ones(NonZeroCPI(t),1)*(theta(17)*a(:,:,jlead1(j)+1)*M(:,:,jlead1(j)+1)+...
        (1-theta(17))*a(:,:,jlead0(j)+1)*M(:,:,jlead0(j)+1))*H;
    Fny=ones(NonZeroNGDP(t),1)*(theta(17)*(a(:,:,jlead1(j)+1)+b(:,:,jlead1(j)+1))*M(:,:,jlead1(j)+1)*H+...
        ((1-theta(17))*(a(:,:,jlead0(j)+1)+b(:,:,jlead0(j)+1))*M(:,:,jlead0(j)+1)*H)-b(:,:,j)*H);
    
    DD=[e1;
        (1-theta(12))*theta(13)*a(:,:,j)+(1-theta(12))*theta(14)*b(:,:,j);
        a(:,:,j);
        b(:,:,j);
        Fpi;
        Fny;];
    
    ZZ=[Z(1:4,t);
        Z(4:3+NonZeroCPI(t),t)-0.25*CPImean;
        Z(maxsurv+4:maxsurv+3+NonZeroNGDP(t),t)-0.25*NGDPmean;];
    
    %         RRR=zeros(size(ZZ,1),size(ZZ,1));
    RRR=[0.0001,zeros(1,NonZeroCPI(t)+2+NonZeroNGDP(t));
        zeros(1,NonZeroCPI(t)+3+NonZeroNGDP(t));
        zeros(2,1),eye(2)*0.01,zeros(2,NonZeroCPI(t)+NonZeroNGDP(t));
        zeros(NonZeroCPI(t),3),eye(NonZeroCPI(t))*((Fpi(1,:)*SigJ(:,:,j)*Fpi(1,:)')^0.5),zeros(NonZeroCPI(t),NonZeroNGDP(t));
        zeros(NonZeroNGDP(t),NonZeroCPI(t)+3),eye(NonZeroNGDP(t))*((Fny(1,:)*SigJ(:,:,j)*Fny(1,:)')^0.5);];
    
    
    Ztilde = ZZ-DD*M(:,:,j)*Xhat(:,t);
    
      Ztilde(2)=Ztilde(2)-theta(12)*Rlag;Rlag=Z(2,t);
        Ztilde(3)=Ztilde(3)-Gr(2,1)*Rlag;
        Ztilde(4)=Ztilde(4)-Gr(1,1)*Rlag;
    PP=M(:,:,j)*pp*M(:,:,j)'+N(:,:,j)*N(:,:,j)';
    OMEGA=DD*PP*DD'+RRR*RRR';
    K=PP*DD'/OMEGA;
    Xhat(:,t+1)= M(:,:,j)*Xhat(:,t)+K*Ztilde;
    pp=PP-K*OMEGA*K';
    
    PP0(:,:,tt)=pp;
    PP1(:,:,tt)=PP;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct what is needed for last step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for tt=1:T;
    s=ST(tt:tt+bindim-1)';
    j=binvec2dec(s)+1;

    Xtilde(:,tt)=M(:,:,j)*Xhat(:,tt);
    
end

Xsm=Xhat*0;
Xsm(:,T)=Xhat(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%backward recursion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tt=T:-1:1
      s=ST(tt:tt+bindim-1)';
    j=binvec2dec(s)+1;
    J=PP0(:,:,tt)*M(:,:,j)'*inv(PP1(:,:,tt)+eye(dimX)*1e-6);
    Xsm(:,tt)=Xhat(:,tt)+J*(Xsm(:,tt+1)-Xtilde(:,tt+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add up and spit out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=Xsm+Xplus;
X=X(:,2:end);