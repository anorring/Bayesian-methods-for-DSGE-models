% Set up and estimate SVAR(1)
clc
clear all
close all

%--------------------------------------------------------------------------
% Set control parameters
%--------------------------------------------------------------------------

%Number of draws
S=1e5;%Number of draws in MCMC
J=300;%Number of draws from MCMC used to simulate posterior functions of theta
USdata=1; %Set to 1 to use US data on inflation and fed funds rate, 0 to use artifical data
epseye=5e-1; %scale up or down to tune acceptance ratio
adaptive=1; %set to 1 to use adaptive proposal density


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Simulate State Space system
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if USdata==0
    
    T=1000;
    c_11=1;
    c_21=-0.5;
    c_22=1;
    
    
    a_11=0.9;
    a_12=-0.2;
    a_21=0.2;
    a_22=0.6;
    
    Ctrue=[c_11,0;
        c_21,c_22;];
    
    Atrue=[a_11,a_12;
        a_21,a_22;];
    
    A=Atrue;
    
    %%
    C=Ctrue;
    U=randn(2,T);
    X=zeros(2,T+1);
    Z=zeros(length(A(:,1)),T);
    for t=1:T;
        X(:,t+1)=A*X(:,t)+C*U(:,t);
        Z(:,t)=X(:,t+1);
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load and detrend US data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load inf-ff
    ff=ff';
    ff(1,:)=detrend(ff(1,:));
    ff(2,:)=detrend(ff(2,:));
    Z=ff;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting values from OLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=length(Z);
Y=Z(:,2:end);
X=Z(:,1:end-1);
beta=Y*X'/(X*X');
Omega_ols= (1/(T-1))*(Y-beta*X)*(Y-beta*X)';
C_ols=chol(Omega_ols)';
A_ols=beta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate posterior using Metroplis Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta=[C_ols(1,1),C_ols(2,1),C_ols(2,2),A_ols(1,1),A_ols(1,2),A_ols(2,1),A_ols(2,2);]';%Starting value for parameter vector
LB=[0,-100,0,-100*ones(1,4);]';
UB=ones(7,1)*100;
x=theta;


%--------------------------------------------------------------------------
% Initializes the MH algorithm
%--------------------------------------------------------------------------

%Initializes the proposal variance.
%Set up so that the first candidate draw is always accepted
lpostdraw = -9e+200;
bdraw=x;
vscale=diag(abs(theta))*1d-3+1e-3*eye(length(x));

thetaMCMC=zeros(length(x),S);%Store all draws in thetaMCMC
%Matrices that keep track of switches and drwas outside LB and UB
OutsideProp=zeros(S,1);
SwitchesProp=zeros(S,1);
%Number of draws outside parameter boundaries
q=0;
%Number of switches (acceptances)
pswitch=0;
%Iteration counter
iter=0;

%--------------------------------------------------------------------------
% MH algorithm starts here
%--------------------------------------------------------------------------

tic
for iter=1:S
    iter;
    % Draw from proposal density Theta*_{t+1} ~ N(Theta_{t},vscale)
    bcan = bdraw + norm_rnd(vscale);
    
    if min(bcan > LB)==1
        if min(bcan < UB)==1
            lpostcan = LLSVAR(bcan,Z);
            laccprob = lpostcan-lpostdraw;
        else
            laccprob=-9e+200;
            q=q+1;
        end
    else
        laccprob=-9e+200;
        q=q+1;
    end
    
    %Accept candidate draw with log prob = laccprob, else keep old draw
    if log(rand)<laccprob
        lpostdraw=lpostcan;
        bdraw=bcan;
        pswitch=pswitch+1;
    end
    
    thetaMCMC(:,iter)=bdraw;
    
    OutsideProp(iter)=q/iter;
    SwitchesProp(iter)=pswitch/iter;
    if adaptive==1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %use for adaptive M-H
    if mod(iter,100)==0;
        vscale=1d-3*cov(thetaMCMC(:,iter)');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
toc

disp(['iter: ',num2str(iter)]);
disp(['acceptance rate: ',num2str(SwitchesProp(iter))]);

figure
thetaMCMC=thetaMCMC(:,2:end);
for jj=1:7;
    subplot(3,3,jj);
    hist(thetaMCMC(jj,:),50);
end

convcheck(thetaMCMC);
plotpost(thetaMCMC,0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find and plot probability intervals for IRFS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
periods=100;
ra = max(size(thetaMCMC));
FF = ceil(ra.*rand(J,1));

IMPR1=zeros(2,periods,J);
IMPR2=zeros(2,periods,J);
VARDECOMP=zeros(2,J);
count=0;
for j=1:J;
    theta=thetaMCMC(:,FF(j,1));
    c_11=theta(1);
    c_21=theta(2);
    c_22=theta(3);
    
    a_11=theta(4);
    a_12=theta(5);
    a_21=theta(6);
    a_22=theta(7);
    
    C=[c_11,0;
        c_21,c_22;];
    
    A=[a_11,a_12;
        a_21,a_22;];
    
    IR1=zeros(2,periods);
    IR2=zeros(2,periods);
    for jj=0:periods-1
        IR1(:,jj+1)= (A^jj)*C(:,1);
        IR2(:,jj+1)= (A^jj)*C(:,2);
    end
    
    count = count +(min(IR2(2,:))<0); %add count if inf on inf is < 0
    
    IMPR1(:,:,j)=IR1;
    IMPR2(:,:,j)=IR2;
    
    %     Variance decomp
    sigx=dlyap(A,C*C');
    sigx1=dlyap(A,C(:,1)*C(:,1)');
    VARDECOMP(:,j)=diag(sigx1)./diag(sigx);
    %
end;



if USdata==0;
    A=Atrue;
    C=Ctrue;
    IR1true=zeros(2,periods);
    IR2true=zeros(2,periods);
    for jj=0:periods-1
        IR1true(:,jj+1)= (A^jj)*C(:,1);
        IR2true(:,jj+1)= (A^jj)*C(:,2);
    end
else
    A=A_ols;
    C=C_ols;
    IR1ols=zeros(2,periods);
    IR2ols=zeros(2,periods);
    for jj=0:periods-1
        IR1ols(:,jj+1)= (A^jj)*C(:,1);
        IR2ols(:,jj+1)= (A^jj)*C(:,2);
    end
end

VARDECOMPsort=sort(VARDECOMP,2);
ImpSort1=sort(IMPR1,3);
ImpSort2=sort(IMPR2,3);

upper=ceil(J*.95);
med=ceil(J*.5);
lower=ceil(J*.05);

figure(5)

subplot(2,2,1);
xlabel('FFR on FFR')
hold on;
plot(ImpSort1(1,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort1(1,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort1(1,:,lower),':','color','black','LineWidth',3);
hold on;
if USdata==0;
    plot(IR1true(1,:),'-','color','red','LineWidth',3);
end

subplot(2,2,2);
xlabel('Inf on FFR')
hold on;
plot(ImpSort1(2,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort1(2,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort1(2,:,lower),':','color','black','LineWidth',3);
hold on;
if USdata==0;
    plot(IR1true(2,:),'-','color','red','LineWidth',3);
end

subplot(2,2,3);
xlabel('FFR on Inf')
hold on;
plot(ImpSort2(1,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort2(1,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort2(1,:,lower),':','color','black','LineWidth',3);
hold on;
if USdata==0;
    plot(IR2true(1,:),'-','color','red','LineWidth',3);
end


subplot(2,2,4);
xlabel('Inf on inf')
hold on;
plot(ImpSort2(2,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort2(2,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort2(2,:,lower),':','color','black','LineWidth',3);
hold on;
if USdata==0;
    plot(IR2true(2,:),'-','color','red','LineWidth',3);
    legend('97.5 %','median','2.5%','True')
else
    legend('97.5 %','median','2.5%')
end

%plot unsorted IRFs
figure(6)
subplot(2,2,1);
plot(reshape(IMPR1(1,:,:),periods,J));
subplot(2,2,2);
plot(reshape(IMPR1(2,:,:),periods,J));
subplot(2,2,3);
plot(reshape(IMPR2(1,:,:),periods,J));
subplot(2,2,4);
plot(reshape(IMPR2(2,:,:),periods,J));

figure(7)
subplot(2,2,1)
hist(VARDECOMP(1,:),25);
subplot(2,2,2);
hist(VARDECOMP(2,:),25);
subplot(2,2,3)
hist(1-VARDECOMP(1,:),25);
subplot(2,2,4);
hist(1-VARDECOMP(2,:),25);

display('Prob(Inf on Inf > 0)')
count/J
