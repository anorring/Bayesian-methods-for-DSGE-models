%counterfactual MDB
clear all;
close all
clc;
format short;
warning off all;
maxsurv=50;
[Z,CPImean,NGDPmean,NonZeroCPI,NonZeroNGDP]=CleanData(maxsurv);
sampleDates=[1981.5:.25:2010.5];

load('lastMCMC');
load('STMCMC');
J=10; % Number of draws; 100 is usually more than enough
ra = max(size(lastMCMC));
f = ceil(ra.*rand(J,1));
T=size(Z,2);
THETAmean=mean(lastMCMC,2);
sdprod= ((1-THETAmean(17))*(THETAmean(3)^2) + THETAmean(17)*((THETAmean(3)*THETAmean(18))^2))^0.5;

%Define hyper parameters etc
kbar=8;
tol=1e-4;
bindim=5;% Number of periods that matter (Markov order, if you like) and dimension of binary identifier of time varying matrices
binmax=2^(bindim); % Number of "regimes"
dimX=kbar*2;

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
Xdist=zeros(dimX,T,J);
Shockdist=zeros(3,T,J);
Xcf=zeros(dimX,T,J);
Zcf=zeros(4,T,J);
for jj=1:J;
    jj
    theta=lastMCMC(:,f(jj,1));
    ST=STMCMC(:,f(jj,1));
    
    
    [P,p,M,N,K,D,L,R,Rj,RRj,SigJ,M0,N0,a0,b0,d0,a,b,d,dimx,dimX,dimu,dimuj,e1,e2,H]= SVLorenz(theta,kbar,tol,binmax);
    [P,p,M,K,N,L,D,R,Rj,RRj,a,b,d,SigJ,EE]  = MOAFLorenz(P,p,M,K,D,N,L,R,Rj,RRj,a,b,d,e1,e2,H,tol,binmax,dimx,dimX,dimu,dimuj,jlag,jlead1,jlead0,SigJ,theta);
    
    [Xss]=MBDsimsmooth(M,N,a,b,e1,bindim,dimX,jlead1,jlead0,SigJ,theta,Z,ST,T,CPImean,NGDPmean,NonZeroCPI,NonZeroNGDP,maxsurv,H);
    Xdist(:,:,jj)=Xss;
    
    Shockdist(1,2:end,jj)=inv(theta(3))*(Z(1,2:end)-M(1,1,1)*Z(1,1:end-1));
    Shockdist(2,2:end,jj)=Xss(2,2:end)-M(2,2,1)*Xss(2,1:end-1);
    Shockdist(3,2:end,jj)=Z(2,2:end)-[theta(13) theta(14)]*Z(3:4,2:end) - theta(12)*Z(2,1:end-1);
    
    thetacf=theta;
    thetacf(17)=0;
    thetacf(18)=1;
    thetacf(3)=((1-theta(17))*(theta(3)^2) + theta(17)*((theta(3)*theta(18))^2))^0.5;
    [P,p,M,N,K,D,L,R,Rj,RRj,SigJ,M0,N0,a0,b0,d0,a,b,d,dimx,dimX,dimu,dimuj,e1,e2,H]= SVLorenz(thetacf,kbar,tol,binmax);
    [P,p,M,K,N,L,D,R,Rj,RRj,a,b,d,SigJ,EE]  = MOAFLorenz(P,p,M,K,D,N,L,R,Rj,RRj,a,b,d,e1,e2,H,tol,binmax,dimx,dimX,dimu,dimuj,jlag,jlead1,jlead0,SigJ,thetacf);
    ST=ST*0;
    
    stick=theta(15); %Calvo parameter
    beta=theta(16); %discount rate
    delta=theta(10); %labour supply curvature
    lambda=(1-stick)*(1-stick*beta)/stick;
    Rlag=theta(12)*Z(2,1);
    Gr=[-1;-lambda+lambda*delta;]*theta(12)/(1-theta(12));
    Gu=[-1;-lambda+lambda*delta;]/(1-theta(12));
     
    if EE==1;
        for tt=2:T;
            s=ST(tt-1:tt-1+bindim-2)';
            j=binvec2dec(s)+1;
            Xcf(:,tt,jj)=M(:,:,j)*Xcf(:,tt-1,jj) + N(:,1:3,j)*Shockdist(1:3,tt,jj);
            Zcf(:,tt,jj)=D(1:4,:,j)*Xcf(:,tt,jj);
            Zcf(2,tt,jj)=Zcf(2,tt,jj)+theta(12)*Rlag;Rlag=Zcf(2,tt,jj);
            Zcf(3,tt,jj)=Zcf(3,tt,jj)+Gr(2,1)*Rlag;
            Zcf(4,tt,jj)=Zcf(4,tt,jj)+Gr(1,1)*Rlag;
            
        end
    else       
                
                Xcf(:,:,jj)=NaN;
                Zcf(:,:,jj)=NaN;
       
    end
end

lower=0.05;
upper=0.95;
ndraws=J;
low=ceil(ndraws*lower);
median=ceil(ndraws*.5);
upp=ceil(ndraws*upper);

Xsort=sort(Xdist,3);
q=ceil(dimX^.5);
MED=reshape(Xsort(:,:,median),dimX,T);
L=reshape(Xsort(:,:,low),dimX,T);
U=reshape(Xsort(:,:,upp),dimX,T);
figure

for jj=1:dimX
    subplot(q,q,jj);
    plot(L(jj,:),':');
    hold on
    plot(MED(jj,:),'k');
    hold on
    plot(U(jj,:),':');
    legend('lower','median','upper');
end

figure(2)

Shocksort=sort(Shockdist,3);
MED=reshape(Shocksort(:,:,median),3,T);
L=reshape(Shocksort(:,:,low),3,T);
U=reshape(Shocksort(:,:,upp),3,T);
figure

for jj=1:3
    subplot(3,1,jj);
    plot(L(jj,:),':');
    hold on
    plot(MED(jj,:),'k');
    hold on
    plot(U(jj,:),':');
    legend('lower','median','upper');
end


figure(3)

Zcfsort=sort(Zcf,3);
MED=reshape(Zcfsort(:,:,median),4,T);
L=reshape(Zcfsort(:,:,low),4,T);
U=reshape(Zcfsort(:,:,upp),4,T);
figure

for jj=1:4
    subplot(2,2,jj);
    plot(Z(jj,:),'r');
    hold on
    plot(L(jj,:),':');
    hold on
    plot(MED(jj,:),'k');
    hold on
    plot(U(jj,:),':');
    legend('actual','lower','median','upper');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot fit at mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
load('xmax')
theta=xmax;
[P,p,M,N,K,D,L,R,Rj,RRj,SigJ,M0,N0,a0,b0,d0,a,b,d,dimx,dimX,dimu,dimuj,e1,e2,H]= SVLorenz(theta,kbar,tol,binmax);
[P,p,M,K,N,L,D,R,Rj,RRj,a,b,d,SigJ,EE]  = MOAFLorenz(P,p,M,K,D,N,L,R,Rj,RRj,a,b,d,e1,e2,H,tol,binmax,dimx,dimX,dimu,dimuj,jlag,jlead1,jlead0,SigJ,theta);

MBD1sidedFit(M,N,a,b,e1,bindim,dimX,jlead1,jlead0,SigJ,theta,Z,ST,T,CPImean,NGDPmean,NonZeroCPI,NonZeroNGDP,maxsurv,H)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot probability weighted innovations to prod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
cutoff=0.9;
load('xmax')
ProdInnov=Z(1,2:end)-xmax(1)*Z(1,1:end-1);
figure
% plot(ProdInnov);
plot(ProdInnov.*(ProdInnov>0));
hold on
plot(abs(ProdInnov.*(ProdInnov<0)),'g');
hold on
% plot(ProdInnov.*mean(STMCMC(7:end,:),2)','k');
figure
plot(ProdInnov.*((mean(STMCMC(7:end,:),2)')>cutoff),'k')

pos=sum((ProdInnov.*((mean(STMCMC(7:end,:),2)')>cutoff)>0));
neg=sum((ProdInnov.*((mean(STMCMC(7:end,:),2)')>cutoff)<0));
MBDsd=sum(abs(ProdInnov).*((mean(STMCMC(7:end,:),2)')>cutoff))/sum(mean(STMCMC(7:end,:),2)'>=cutoff);
sd=mean(abs(ProdInnov));

[pos neg sd MBDsd]



