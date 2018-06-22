function [LL]=MBDLL_Lorenz(M,N,a,b,e1,dimX,SigJ,theta,Z,T,CPImean,NGDPmean,NonZeroCPI,NonZeroNGDP,maxsurv,H)

% try
    pp=M*0;
    for j=1:10;
        pp=M*pp*M'+N*N';
    end
    X=zeros(dimX,T);
stick=theta(16); %Calvo parameter
beta=theta(17); %discount rate
varphi=theta(11); %labour supply curvature
omega=1;
lambda=(1-stick)*(1-stick*beta)/stick;

    LL=0;
    Rlag=theta(13)*Z(2,1);
    Gr=[-1;-lambda+lambda*varphi;]*theta(13)/(1-theta(13));
%     Gu=[-1;-lambda+lambda*varphi;]/(1-theta(13));
    for t=1:T-1;
        
        %%
        
        fpi=a*M;
        fny=(a+b)*M - b;
        Fpi=ones(NonZeroCPI(t),1)*fpi*H;
        Fny=ones(NonZeroNGDP(t),1)*fny*H;
        
        DD=[e1;
            (1-theta(13))*theta(14)*a+(1-theta(14))*theta(15)*b;
            a;
            b;
            Fpi;
            Fny;];
        
        ZZ=[Z(1:4,t);
            Z(4:3+NonZeroCPI(t),t)-0.25*CPImean;
            Z(maxsurv+4:maxsurv+3+NonZeroNGDP(t),t)-0.25*NGDPmean;];
        

        RRR=[0.001,zeros(1,NonZeroCPI(t)+3+NonZeroNGDP(t));
            0,theta(5),zeros(1,NonZeroCPI(t)+2+NonZeroNGDP(t)); 
           
            zeros(2,1),[Gr(2);Gr(1);]*theta(5), eye(2)*0.001,zeros(2,NonZeroCPI(t)+NonZeroNGDP(t));
            zeros(NonZeroCPI(t),4),eye(NonZeroCPI(t))*((fpi(1,:)*SigJ*fpi(1,:)')^0.5),zeros(NonZeroCPI(t),NonZeroNGDP(t));
            zeros(NonZeroNGDP(t),NonZeroCPI(t)+4),eye(NonZeroNGDP(t))*((fny(1,:)*SigJ*fny(1,:)')^0.5);];
         
        
        Ztilde = ZZ-DD*M*X(:,t);
        Ztilde(2)=Ztilde(2)-theta(13)*Rlag;Rlag=Z(2,t);
        Ztilde(3)=Ztilde(3)-Gr(2,1)*Rlag;
        Ztilde(4)=Ztilde(4)-Gr(1,1)*Rlag;
        PP=M*pp*M'+N*N';
        OMEGA=DD*PP*DD'+RRR*RRR';
        
        LL=LL-0.5*(logdet(OMEGA)+Ztilde'/OMEGA*Ztilde);       
        
        K=PP*DD'/OMEGA;
        X(:,t+1)= M*X(:,t)+K*Ztilde;
        pp=PP-K*OMEGA*K';
        
        
        
    end
    
        
    if isnan(LL)==1;
        LL=-1e200;
         display('LL is NaN')
    end
    
    if LL==Inf;
        LL=-1e200;
         display('LL is infinite')
    end
    if isreal(LL)==0;
        LL=-1e200;
        display('LL is imaginary')
    end
    
% catch
%     LL=-1e200;
%     display('Not possible to evaluate LL')
% 
% end

