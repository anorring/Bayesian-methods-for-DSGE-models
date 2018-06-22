clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% Example: man bites dog
rho1=0.87; %persistence of technology
rho2=0.72; %persistence of demand
sigu=.018;  %s.d. of state innov
sigud=0.01; %s.d. "demand" shock
sigur=.014;  %s.d. of m.p. shock
sigaj=.25;  %s.d. of island tech
sigzj1=.4;  %s.d. of private info noise
sigzj2=.4;  %s.d. of private info noise
sigdj=0.21;  %s.d. of island demand
sigmbd=0.5; %s.d. of m-b-d signal
varphi=1.1; %labour supply curvature
delta=1.07; %elasticity of demand
fir=0.08;%Interest inertia
fipi=1.55; %Taylor param;
fiy=0.24; %Taylor rule param
stick=0.74; %Calvo parameter
beta=0.99; %discount rate
omega=0.05;  %unconditional prob of S=1, i.e. of observing pub signal
gamma=6.27;    %s.d. multiplier of u when S=1

theta=[rho1,rho2,sigu,sigud,sigur,sigaj,sigzj1,sigzj2,sigdj,sigmbd,varphi,delta,fir,fipi,fiy,stick,beta,omega,gamma;];
load('xmax');
theta=xmax;
% theta(12)=0;
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

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign starting values to endogenous matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%
[P,p,M,N,K,D,L,R,Rj,RRj,SigJ,M0,N0,a0,b0,d0,a,b,d,dimx,dimX,dimu,dimuj,e1,e2,H]= SVLorenz(theta,kbar,tol,binmax);
toc
%%
% s=3;
% for t=1:20
%     infl(t)=a0*M0^(t-1)*N0(:,s);
%     output(t)=b0*M0^(t-1)*N0(:,s);
%     hier(:,t)=M0^(t-1)*N0(:,s);
%     
% end
% 
% 
% subplot(3,1,1);plot(infl);
% subplot(3,1,2);plot(output);
% 
% subplot(3,1,3);plot(hier(1:5,:)');
%%
tic
[P,p,M,K,N,L,D,R,Rj,RRj,a,b,d,SigJ,EE]  = MOAFLorenz(P,p,M,K,D,N,L,R,Rj,RRj,a,b,d,e1,e2,H,tol,binmax,dimx,dimX,dimu,dimuj,jlag,jlead1,jlead0,SigJ,theta);
toc
%%
%IRF
figure
periods=25;
X1=zeros(dimX,periods);
X0=zeros(dimX,periods);

shock=1;
X0(:,1)=N(:,shock,1);
X1(:,1)=(1/gamma)*N(:,shock,2);
% X1(:,1)=N(:,shock,2);
infl1(1)=a(:,:,2)*X1(:,1);
infl0(1)=a(:,:,1)*X0(:,1);
outp1(1)=b(:,:,2)*X1(:,1);
outp0(1)=b(:,:,1)*X0(:,1);
for t=1:periods;
    s=zeros(1,bindim);
    if t <= bindim;
    s(t)=1;
    end
    s=s(end:-1:1);
    j=binvec2dec(s);
    X1(:,t+1)=M(:,:,j+1)*X1(:,t);
    infl1(t+1)=a(1,:,j+1)*X1(:,t+1);
    outp1(t+1)=b(1,:,j+1)*X1(:,t+1);
    mc1(t)=(1+delta)*(b(1,:,j+1)-e1)*X1(:,t+1);
    X0(:,t+1)=M(:,:,1)*X0(:,t);
    infl0(t+1)=a(1,:,1)*X0(:,t+1);
    outp0(t+1)=b(1,:,1)*X0(:,t+1);
    
end
% 
% shock=2;
% X0(:,1)=N(:,shock,1);
% X1(:,1)=(1/gamma)*N(:,shock,2);
% infl1(1)=a(:,:,2)*X1(:,1);
% infl0(1)=a(:,:,1)*X0(:,1);
% outp1(1)=b(:,:,2)*X1(:,1);
% outp0(1)=b(:,:,1)*X0(:,1);
% for t=1:periods;
%     s=zeros(1,bindim);
%     if t <= bindim;
%     s(t)=1;
%     end
%     s=s(end:-1:1);
%     j=binvec2dec(s);
%     X1(:,t+1)=M(:,:,j+1)*X1(:,t);
%     infl2(t+1)=a(1,:,j+1)*X1(:,t+1);
%     outp2(t+1)=b(1,:,j+1)*X1(:,t+1);
%     mc1(t)=(1+delta)*(b(1,:,j+1)-e1)*X1(:,t+1);
%     X0(:,t+1)=M(:,:,1)*X0(:,t);
%     infl3(t+1)=a(1,:,1)*X0(:,t+1);
%     outp3(t+1)=b(1,:,1)*X0(:,t+1);
%     
% end

subplot(2,1,1);plot(infl1); hold on;plot(infl0,':');
subplot(2,1,2);plot(outp1); hold on;plot(outp0,':');
% subplot(4,1,3);plot(infl2); hold on;plot(infl3,':');
% subplot(4,1,4);plot(outp2); hold on;plot(outp3,':');
% subplot(4,1,3);plot(X0(2:2:12,:)'); 
% subplot(4,1,4);plot(X1(2:2:12,:)');


% %%
% %ARCH
% periods=10;
% MULTinfl=zeros(1,periods);
% MULToutp=zeros(1,periods);
% 
% for t=1:periods;
%     s=zeros(1,bindim);
%     if t <= bindim;
%     s(t)=1;
%     end
%     s=s(end:-1:1);
%     j=binvec2dec(s);
%     if t==1;
%     MULTIinfl(t)=(1/gamma)*a(1,:,j+1)*N(:,shock,j+1);
%     MULTIoutp(t)=(1/gamma)*b(1,:,j+1)*N(:,shock,j+1);
%     else
%     
%     MULTIinfl(t)=a(1,:,j+1)*N(:,shock,j+1);
%     MULTIoutp(t)=b(1,:,j+1)*N(:,shock,j+1);
%        
%     end
%     
% end
% figure
% plot(MULTIinfl); 
% figure
% plot(MULTIoutp);
% 
% %%
% %Dispersion
% periods=bindim;
% DISP=zeros(1,periods);
% 
% 
% for t=1:periods;
%     s=zeros(1,bindim);
%     if t <= bindim;
%     s(t)=1;
%     end
%     s=s(end:-1:1);
%     j=binvec2dec(s);
%     
%     DISP(t)=(a(1,:,jlead0(j)+1)*M(:,:,jlead0(j)+1))*SigJ(:,:,j+1)*(a(1,:,jlead0(j)+1)*M(:,:,jlead0(j)+1))';
%     
%     
% end
% figure
% plot(DISP);
% xlabel('Cross sectional variance of inflation expectations');
% 
% 
% 
% 
