function [LL]=MBDLL(M,N,a,b,e1,bindim,dimX,jlead1,jlead0,SigJ,theta,Z,ST,T,CPImean,NGDPmean,NonZeroCPI,NonZeroNGDP,maxsurv,H)

 try
    pp=M(:,:,2)*0;
    for j=1:10;
        pp=M(:,:,2)*pp*M(:,:,2)'+N(:,:,2)*N(:,:,2)';
    end
    X=zeros(dimX,T);
stick=theta(16); %Calvo parameter
beta=theta(17); %discount rate
varphi=theta(11); %labour supply curvature

lambda=(1-stick)*(1-stick*beta)/stick;

    LL=0;
    Rlag=theta(13)*Z(2,1);
    Gr=[-1;-lambda+lambda*varphi;]*theta(13)/(1-theta(13));
    Gu=[-1;-lambda+lambda*varphi;]/(1-theta(13));
    for t=1:T-1;
        
        %%
        
        s=ST(t:t+bindim-1)';
        j=binvec2dec(s)+1;
        fpi=(theta(18)*a(:,:,jlead1(j)+1)*M(:,:,jlead1(j)+1)+...
            (1-theta(18))*a(:,:,jlead0(j)+1)*M(:,:,jlead0(j)+1));
        Fpi=ones(NonZeroCPI(t),1)*fpi*H;
        fny=theta(18)*(a(:,:,jlead1(j)+1)+b(:,:,jlead1(j)+1))*M(:,:,jlead1(j)+1)+...
            (1-theta(18))*(a(:,:,jlead0(j)+1)+b(:,:,jlead0(j)+1))*M(:,:,jlead0(j)+1);
        Fny=ones(NonZeroNGDP(t),1)*(fny-b(:,:,j))*H;
        
        DD=[e1;
            (1-theta(13))*theta(14)*a(:,:,j)+(1-theta(14))*theta(15)*b(:,:,j);
            a(:,:,j);
            b(:,:,j);
            Fpi;
            Fny;];
        
        ZZ=[Z(1:4,t);
            Z(4:3+NonZeroCPI(t),t)-0.25*CPImean;
            Z(maxsurv+4:maxsurv+3+NonZeroNGDP(t),t)-0.25*NGDPmean;];
        
        %         RRR=zeros(size(ZZ,1),size(ZZ,1));
        RRR=[0.001,zeros(1,NonZeroCPI(t)+3+NonZeroNGDP(t));
            0,theta(5),zeros(1,NonZeroCPI(t)+2+NonZeroNGDP(t)); 
           
            zeros(2,1),[Gr(2);Gr(1);]*theta(5), eye(2)*0.001,zeros(2,NonZeroCPI(t)+NonZeroNGDP(t));
            zeros(NonZeroCPI(t),4),eye(NonZeroCPI(t))*((fpi*SigJ(:,:,j)*fpi')^0.5),zeros(NonZeroCPI(t),NonZeroNGDP(t));
            zeros(NonZeroNGDP(t),NonZeroCPI(t)+4),eye(NonZeroNGDP(t))*((fny*SigJ(:,:,j)*fny')^0.5);];
         
        
        Ztilde = ZZ-DD*M(:,:,j)*X(:,t);
        Ztilde(2)=Ztilde(2)-theta(13)*Rlag;Rlag=Z(2,t);
        Ztilde(3)=Ztilde(3)-Gr(2,1)*Rlag;
        Ztilde(4)=Ztilde(4)-Gr(1,1)*Rlag;
        PP=M(:,:,j)*pp*M(:,:,j)'+N(:,:,j)*N(:,:,j)';
        OMEGA=DD*PP*DD'+RRR*RRR';
        OMEGAinv=eye(length(Ztilde))/OMEGA;
        LL=LL-0.5*(logdet(OMEGA)+Ztilde'*OMEGAinv*Ztilde);   
        
%          if LL==Inf;
%          LL=LL-0.5*(log(det(OMEGA+1e-6*eye(length(OMEGA))))+Ztilde'/(OMEGA+1e-6*eye(length(OMEGA)))*Ztilde);
%          end
        
        K=PP*DD'*OMEGAinv;
        X(:,t+1)= M(:,:,j)*X(:,t)+K*Ztilde;
        pp=PP-K*OMEGA*K';
        
        
        
    end
    
    if LL==Inf;
        LL=-1e200;
    end
    if isreal(LL)==0;
        LL=-1e200;
    end
    
catch
    LL=-1e200;
    display('Not possible to evaluate LL')

end

