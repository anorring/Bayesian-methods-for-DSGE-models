%Program to construct graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all
clc;
format short;
warning off all;
load('Qdoms')
load('xmax')
maxsurv=25;
[Z,CPImean,NGDPmean,NonZeroCPI,NonZeroNGDP]=CleanData(maxsurv);

for r=1:size(Z,1)
    for c=1:size(Z,2)
        
        if Z(r,c)==0;
            Z(r,c)=NaN;
        end
    end
end

sampleDates=[1981.5:.25:2010.5];

figure(1)
subplot(4,1,1);plot(sampleDates,Z(1,:),'linewidth',2);
subplot(4,1,2);plot(sampleDates,Z(2,:),'linewidth',2);
subplot(4,1,3);plot(sampleDates,Z(3,:),'linewidth',2);
subplot(4,1,4);plot(sampleDates,Z(4,:),'linewidth',2);

figure(2)
subplot(2,1,1);plot(sampleDates,Z(5:maxsurv+4,:)-0.25*CPImean);
hold on;
subplot(2,1,1);plot(sampleDates,Z(3,:),'linewidth',2,'color','black','linestyle','--');
subplot(2,1,2);plot(sampleDates(2:end),Z(5+maxsurv:end,2:end)-0.25*NGDPmean);
hold on;
subplot(2,1,2);plot(sampleDates(2:end),Z(3,2:end)+(Z(4,2:end)-Z(4,1:end-1)),'linewidth',2,'color','black','linestyle','--');

%%
%Plot posterior estiamte of s(T)
figure(3)
load('STMCMC');
% STMCMC=STMCMC(:,1500:end);
subplot(2,1,1);
plot(sampleDates,mean(STMCMC(6:end,:),2)*NaN,'linewidth',3,'color','blue');
Rdates=[1981.5 1982.75; 1990.5 1991; 2001 2001.75; 2007.75:2009.25; ];
ycord=ylim;
recessionbars(Rdates,ycord)
hold on
plot(sampleDates,mean(STMCMC(6:end,:),2),'linewidth',3,'color','blue');
hold on
plot(sampleDates,Qdoms,'linewidth',1,'color','red');


%%
% Plot posterior of impact of unit prod shock

load('lastMCMC');
% lastMCMC=lastMCMC(:,1500:end);
J=100; % Number of draws; 100 is usually more than enough
ra = max(size(lastMCMC));
f = ceil(ra.*rand(J,1));
T=size(STMCMC,1);
THETAmean=mean(lastMCMC,2);
sdprod= ((1-THETAmean(17))*(THETAmean(3)^2) + THETAmean(17)*((THETAmean(3)*THETAmean(18))^2))^0.5;

  %Define hyper parameters etc
    kbar=8;
    tol=1e-4;
    bindim=5;% Number of periods that matter (Markov order, if you like) and dimension of binary identifier of time varying matrices
    binmax=2^(bindim); % Number of "regimes"
    dimS=2;
    for j=1:bindim;
        binbase(1,bindim-j+1)=2^(j-1);
    end
    
    jlag=zeros(binmax,1);
    jlead1=zeros(binmax,1);
    jlead0=zeros(binmax,1);
    for j=1:binmax;
        zvec=dec2binvec((j-1),bindim);
        jlag(j)=binbase*[0 zvec(1,1:end-1) ]';
        jlead1(j)=binbase*[zvec(1,2:end) 1 ]';
        jlead0(j)=binbase*[zvec(1,2:end) 0 ]';
    end
    periods=8;
PRODIMP=zeros(J,T-bindim+1);
ARCH=zeros(J,periods);
IMP=zeros(J,100);
for jj=1:J;
    jj
    theta=lastMCMC(:,f(jj,1));
    s=STMCMC(:,f(jj,1));
    
  
    [P,p,M,N,K,D,L,R,Rj,RRj,SigJ,M0,N0,a0,b0,d0,a,b,d,dimx,dimX,dimu,dimuj,e1,e2,H]= SVLorenz(theta,kbar,tol,binmax);
    [P,p,M,K,N,L,D,R,Rj,RRj,a,b,d,SigJ,EE]  = MOAFLorenz(P,p,M,K,D,N,L,R,Rj,RRj,a,b,d,e1,e2,H,tol,binmax,dimx,dimX,dimu,dimuj,jlag,jlead1,jlead0,SigJ,theta);
 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Time series of prod impact
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t=bindim:T
        st=s(t-bindim+1:t);
        j=binvec2dec(st');
        PRODIMP(jj,t)=sdprod*b(:,:,j+1)*N(:,1,j+1)/N(1,1,j+1);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% ARCH and cross sectional dispersion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  SVEC=[ones(1, bindim-1)  1 zeros(1,periods-1) ];

for t=1:periods;
    st=SVEC(1,t:t+bindim-1);
    j=binvec2dec(st);
    ARCH(jj,t)=sdprod*b(:,:,j+1)*N(:,1,j+1)/N(1,1,j+1);
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%Conditional Impact
%%%%%%%%%%%%%%%%%%%%%%%%%   
incr=0.05;
for t=1:100;
    
    IMP0=sdprod*b(:,:,1)*N(:,1,1)/N(1,1,1);
    IMP1=sdprod*b(:,:,2)*N(:,1,2)/N(1,1,2);
    sig=sdprod*t*incr;
    sigplot(t)=sig;
    g=theta(17);
    weight1=(g/((2*pi*(theta(18)*theta(3))^2)^.5))*exp(-(sig^2)/(2*((theta(18)*theta(3))^2)));
    weight2=((1-g)/((2*pi*(theta(3))^2)^.5))*exp(-(sig^2)/(2*(theta(3)^2)));
    weight(t)=weight1/(weight1+weight2);
    IMP(jj,t)=weight(t)*IMP1 + (1-weight(t))*IMP0;
end  


%%%%%%%%%%%%%%%%%%%%%%%%%
%Impulse
%%%%%%%%%%%%%%%%%%%%%%%%%
periods=15;
X1=zeros(dimX,periods);
X0=zeros(dimX,periods);
shock=1;
X0(:,1)=sdprod*N(:,shock,1)/N(1,shock,1);
X1(:,1)=sdprod*N(:,shock,2)/N(1,shock,2);
infl1(1)=a(:,:,2)*X1(:,1);
infl0(1)=a(:,:,1)*X0(:,1);
outp1(1)=b(:,:,2)*X1(:,1);
outp0(1)=b(:,:,1)*X0(:,1);

 SVEC=[ones(1, bindim-1) 1 0 0 0 zeros(1,periods-3) ];
 
 
for t=1:periods;
    st=SVEC(1,t:t+bindim-1);
    j=binvec2dec(st);
    X1(:,t+1)=M(:,:,j+1)*X1(:,t);
    infl1(t+1)=a(1,:,j+1)*X1(:,t+1);
    outp1(t+1)=b(1,:,j+1)*X1(:,t+1);
    X0(:,t+1)=M(:,:,1)*X0(:,t);
    infl0(t+1)=a(1,:,1)*X0(:,t+1);
    outp0(t+1)=b(1,:,1)*X0(:,t+1);
    
end

INFL0(jj,:)=infl0;
INFL1(jj,:)=infl1;
OUTP0(jj,:)=outp0;
OUTP1(jj,:)=outp1;

%Lorenzoni (public signal always there)
% thetacf=xmax;
    thetacf(17)=0;
    thetacf(18)=1;
    thetacf(3)=((1-theta(17))*(theta(3)^2) + theta(17)*((theta(3)*theta(18))^2))^0.5;
    [P,p,M,N,K,D,L,R,Rj,RRj,SigJ,M0,N0,a0,b0,d0,a,b,d,dimx,dimX,dimu,dimuj,e1,e2,H]= SVLorenz(thetacf,kbar,tol,binmax);
    [P,p,M,K,N,L,D,R,Rj,RRj,a,b,d,SigJ,EE]  = MOAFLorenz(P,p,M,K,D,N,L,R,Rj,RRj,a,b,d,e1,e2,H,tol,binmax,dimx,dimX,dimu,dimuj,jlag,jlead1,jlead0,SigJ,thetacf);
   
periods=15;
X1=zeros(dimX,periods);
X0=zeros(dimX,periods);
shock=1;
X0(:,1)=sdprod*N(:,shock,1)/N(1,shock,1);
X1(:,1)=sdprod*N(:,shock,2)/N(1,shock,2);
infl1(1)=a(:,:,2)*X1(:,1);
infl0(1)=a(:,:,1)*X0(:,1);
outp1(1)=b(:,:,2)*X1(:,1);
outp0(1)=b(:,:,1)*X0(:,1);

 SVEC=[ones(1, bindim-1) 1 0 0 0 zeros(1,periods-3) ];
 
 
for t=1:periods;
    st=SVEC(1,t:t+bindim-1);
    j=binvec2dec(st);
    X1(:,t+1)=M(:,:,j+1)*X1(:,t);
    infl1(t+1)=a(1,:,j+1)*X1(:,t+1);
    outp1(t+1)=b(1,:,j+1)*X1(:,t+1);
    X0(:,t+1)=M(:,:,1)*X0(:,t);
    infl0(t+1)=a(1,:,1)*X0(:,t+1);
    outp0(t+1)=b(1,:,1)*X0(:,t+1);
    
end
LORINFL0(jj,:)=infl0;
LORINFL1(jj,:)=infl1;
LOROUTP0(jj,:)=outp0;
LOROUTP1(jj,:)=outp1;

%Lorenzoni (random arrival of signal)
% thetacf=xmax;
%     thetacf(17)=0;
    thetacf(18)=1;
    thetacf(3)=((1-theta(17))*(theta(3)^2) + theta(17)*((theta(3)*theta(18))^2))^0.5;
    [P,p,M,N,K,D,L,R,Rj,RRj,SigJ,M0,N0,a0,b0,d0,a,b,d,dimx,dimX,dimu,dimuj,e1,e2,H]= SVLorenz(thetacf,kbar,tol,binmax);
    [P,p,M,K,N,L,D,R,Rj,RRj,a,b,d,SigJ,EE]  = MOAFLorenz(P,p,M,K,D,N,L,R,Rj,RRj,a,b,d,e1,e2,H,tol,binmax,dimx,dimX,dimu,dimuj,jlag,jlead1,jlead0,SigJ,thetacf);
   
periods=15;
X1=zeros(dimX,periods);
X0=zeros(dimX,periods);
shock=1;
X0(:,1)=sdprod*N(:,shock,1)/N(1,shock,1);
X1(:,1)=sdprod*N(:,shock,2)/N(1,shock,2);
infl1(1)=a(:,:,2)*X1(:,1);
infl0(1)=a(:,:,1)*X0(:,1);
outp1(1)=b(:,:,2)*X1(:,1);
outp0(1)=b(:,:,1)*X0(:,1);

 SVEC=[ones(1, bindim-1) 1 0 0 0 zeros(1,periods-3) ];
 
 
for t=1:periods;
    st=SVEC(1,t:t+bindim-1);
    j=binvec2dec(st);
    X1(:,t+1)=M(:,:,j+1)*X1(:,t);
    infl1(t+1)=a(1,:,j+1)*X1(:,t+1);
    outp1(t+1)=b(1,:,j+1)*X1(:,t+1);
    X0(:,t+1)=M(:,:,1)*X0(:,t);
    infl0(t+1)=a(1,:,1)*X0(:,t+1);
    outp0(t+1)=b(1,:,1)*X0(:,t+1);
    
end
RANDLORINFL0(jj,:)=infl0;
RANDLORINFL1(jj,:)=infl1;
RANDLOROUTP0(jj,:)=outp0;
RANDLOROUTP1(jj,:)=outp1;

end
% plot(sigplot./sdprod,IMP)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct posterior probability intervals for theta
% Used in Table 1 of paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lower=ceil(J*.025);
median=ceil(J*.5);
upper=ceil(J*.975);

%
PRODIMPsort=sort(PRODIMP,1);
PRODIMPlow=PRODIMPsort(lower,:);
PRODIMPmed=PRODIMPsort(median,:);
PRODIMPupp=PRODIMPsort(upper,:);


INFL0sort=sort(INFL0,1);
INFL1sort=sort(INFL1,1);
OUTP0sort=sort(OUTP0,1);
OUTP1sort=sort(OUTP1,1);

LORINFL0sort=sort(LORINFL0,1);
LORINFL1sort=sort(LORINFL1,1);
LOROUTP0sort=sort(LOROUTP0,1);
LOROUTP1sort=sort(LOROUTP1,1);

RANDLORINFL0sort=sort(RANDLORINFL0,1);
RANDLORINFL1sort=sort(RANDLORINFL1,1);
RANDLOROUTP0sort=sort(RANDLOROUTP0,1);
RANDLOROUTP1sort=sort(RANDLOROUTP1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(2,1,2);
plot(sampleDates,PRODIMPmed(1,6:end)*NaN,'linestyle','-','linewidth',3,'color','blue');

Rdates=[1981.5 1982.75; 1990.5 1991; 2001 2001.75; 2007.75:2009.25; ];
v=AXIS;
ycord=[v(3),v(4)];
recessionbars(Rdates,ycord)
hold on

plot(sampleDates,PRODIMPlow(1,6:end),'linestyle',':','linewidth',2,'color','blue');
hold on
plot(sampleDates,PRODIMPmed(1,6:end),'linestyle','-','linewidth',3,'color','blue');
hold on
plot(sampleDates,PRODIMPupp(1,6:end),'linestyle',':','linewidth',2,'color','blue');


figure(4)
ARCHsort=sort(ARCH,1);
plot(ARCHsort(lower,1:6),'linestyle',':','linewidth',2,'color','blue');
hold on
plot(ARCHsort(median,1:6),'linestyle','-','linewidth',3,'color','blue');
hold on
plot(ARCHsort(upper,1:6),'linestyle',':','linewidth',2,'color','blue');

figure(5)
IMPsort=sort(IMP,1);
plot([incr:incr:incr*100],IMPsort(lower,:),'linestyle',':','linewidth',2,'color','blue');
hold on
plot([incr:incr:incr*100],IMPsort(median,:),'linestyle','-','linewidth',3,'color','blue');
hold on
plot([incr:incr:incr*100],IMPsort(upper,:),'linestyle',':','linewidth',2,'color','blue');

figure(6)
subplot(2,1,1);
plot(INFL0sort(lower,:),'linestyle',':','linewidth',2,'color','black');
hold on
plot(INFL0sort(median,:),'linestyle','-','linewidth',3,'color','black');
hold on
plot(INFL0sort(upper,:),'linestyle',':','linewidth',2,'color','black');
hold on
plot(INFL1sort(lower,:),'linestyle',':','linewidth',2,'color','blue');
hold on
plot(INFL1sort(median,:),'linestyle','-','linewidth',3,'color','blue');
hold on
plot(INFL1sort(upper,:),'linestyle',':','linewidth',2,'color','blue');

subplot(2,1,2);
plot(OUTP0sort(lower,:),'linestyle',':','linewidth',2,'color','black');
hold on
plot(OUTP0sort(median,:),'linestyle','-','linewidth',3,'color','black');
hold on
plot(OUTP0sort(upper,:),'linestyle',':','linewidth',2,'color','black');
hold on
plot(OUTP1sort(lower,:),'linestyle',':','linewidth',2,'color','blue');
hold on
plot(OUTP1sort(median,:),'linestyle','-','linewidth',3,'color','blue');
hold on
plot(OUTP1sort(upper,:),'linestyle',':','linewidth',2,'color','blue');

%Lorenzoni
subplot(2,1,1);
hold on
plot(LORINFL1sort(lower,:),'linestyle',':','linewidth',2,'color','green');
hold on
plot(LORINFL1sort(median,:),'linestyle','-','linewidth',3,'color','green');
hold on
plot(LORINFL1sort(upper,:),'linestyle',':','linewidth',2,'color','green');


subplot(2,1,2);
hold on
plot(LOROUTP1sort(lower,:),'linestyle',':','linewidth',2,'color','green');
hold on
plot(LOROUTP1sort(median,:),'linestyle','-','linewidth',3,'color','green');
hold on
plot(LOROUTP1sort(upper,:),'linestyle',':','linewidth',2,'color','green');




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %IRF of cross-sectional dispersion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta=mean(lastMCMC,2);


    [P,p,M,N,K,D,L,R,Rj,RRj,SigJ,M0,N0,a0,b0,d0,a,b,d,dimx,dimX,dimu,dimuj,e1,e2,H]= SVLorenz(theta,kbar,tol,binmax);
    [P,p,M,K,N,L,D,R,Rj,RRj,a,b,d,SigJ,EE]  = MOAFLorenz(P,p,M,K,D,N,L,R,Rj,RRj,a,b,d,e1,e2,H,tol,binmax,dimx,dimX,dimu,dimuj,jlag,jlead1,jlead0,SigJ,theta);
 

periods=20;
SVEC=[zeros(1, bindim) 1 zeros(1,periods-1) ]; %define sequence of Z
 st=SVEC(1,1:bindim);
    j=binvec2dec(st);
    
LM=(eye(dimX)-K(:,:,j+1)*D(:,:,j+1))*M(:,:,j+1);

ISV=dlyap(LM,K(:,:,j+1)*Rj(:,:,j+1)*Rj(:,:,j+1)'*K(:,:,j+1)');
DISP1=[];
DISP1(1,1)= e1*ISV*(e1)';
SIG=ISV;


for t=1:periods;
    st=SVEC(1,t:t+bindim-1);
    j=binvec2dec(st);
    
    
        LM=(eye(dimX)-K(:,:,j+1)*D(:,:,j+1))*M(:,:,j+1);
        KRRK=K(:,:,j+1)*Rj(:,:,j+1)*Rj(:,:,j+1)'*K(:,:,j+1)';
        SIG=LM*SIG*LM'+ KRRK;

DISP1(1,t+1)= e1*SIG*(e1)';
end  
 %%
figure(7)
plot(DISP1)
xx=[];
%normal plotter
fDISP1=zeros(200,10);
axmin=max(DISP1.^.5)*.01;
for jj=1:10;
for j=1:200;
 sig=DISP1(1,jj); 
   x=-axmin+j*axmin*0.01;
 
    fy(j,1)=(1/((sig*(2*pi)^.5)))*exp(-(x^2)/(2*sig^2));  
    fDISP1(j,jj)=fy(j,1);
   
end


end
figure(8)
surf([1:10],[-axmin:0.01*axmin:axmin-0.01*axmin],fDISP1);

J=length(lastMCMC);
lower=ceil(J*.025);
median=ceil(J*.5);
upper=ceil(J*.975);
MCMCsort=sort(lastMCMC,2);
thetalow=MCMCsort(:,lower);
thetamed=MCMCsort(:,median);
thetaupper=MCMCsort(:,upper);
