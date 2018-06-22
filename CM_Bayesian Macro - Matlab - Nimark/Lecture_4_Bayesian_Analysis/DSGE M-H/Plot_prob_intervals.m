%Program to plot probability intervals for IRFs etc using MCMC from M-H
%chain.
clc
clear all
close all
global Z
load('Z');
load('bb_');
load('bb_mle');
load('bb_prior');
%  bb_=bb_prior;%Uncomment for prior predicitive analysis
%  bb_=bb_mle;%Uncomment for improper uniform priors
burnin=0.5*length(bb_);
bb_=bb_(:,burnin:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=2000;% NUmber of independent draws from the markov chain
periods=25; %Number of periods in the IRFs
T=length(Z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find and plot probability intervals for IRFS, Variance decompisition and
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ra = max(size(bb_));
ff = ceil(ra.*rand(S,1));
count=0;

IMPR_r=zeros(periods,S);%This is where we store the IRFs for the different draws
IMPR_pi=zeros(periods,S);
IMPR_y=zeros(periods,S);
VARDECOMP=zeros(3,S);%This is where we store the variance decompositions for the different draws
forecast=zeros(3,periods+1,S);
for s=1:S;
    theta=bb_(:,ff(s,1));%Take a random draw from the MCMC
    [A,C,D,R]=DSGE_SS(theta);% Compute the state space matrices for theta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Impulse response functiom
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    IR=zeros(periods,3);
    for ss=0:periods-1
        IR(ss+1,:)= D*(A^ss)*C(1,1);
    end
    IMPR_r(:,s)=IR(:,1);
    IMPR_pi(:,s)=IR(:,2);
    IMPR_y(:,s)=IR(:,3);
    %%%%%%%%%%%%%%%%%%%%%%%%
    %Variance decomp
    %%%%%%%%%%%%%%%%%%%%%%%%
    sigx=theta(6);%Picks out s.d. of prod innovations from theta
    SIG_X=sigx^2*inv(1-A^2);%Unconditional variance of productivity
    RR=R*R';% variance due to moneetary policy, cost push and demand shocks
    SIG_Y= D*SIG_X*D' +RR;% total unconditional
    
    VARDECOMP(:,s)=ones(3,1) - diag(RR)./diag(SIG_Y); % Compute fraction of unconditional variance due to productivity shocks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Count instances of logical statement being true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  IR(1,3) > 0.005;
        count=count+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Make forecasts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [K,P,p,Xtt]=kalman_filter(A,C,D,R,Z);
    X=Xtt(:,end)+randn*(p^0.5);
    Xst=X;
    forecast(:,1,s)=Z(:,end);
    for ss=2:periods+1
        X=A*X+C(1,1)*randn;
    forecast(:,ss,s)=D*X+ R*randn(4,1);
    end
end;

VARDECOMPsort=sort(VARDECOMP,2);
ImpSort_r=sort(IMPR_r,2);
ImpSort_pi=sort(IMPR_pi,2);
ImpSort_y=sort(IMPR_y,2);

figure(1)
subplot(3,1,1);
hold on;
plot(ImpSort_r(:,S*.95),':','color','black','LineWidth',3);
hold on;
plot(ImpSort_r(:,S*.5),'color','black','LineWidth',3);
hold on;
plot(ImpSort_r(:,S*.05),':','color','black','LineWidth',3);
xlabel('Interest rate response to productivity shock')

subplot(3,1,2);
hold on;
plot(ImpSort_pi(:,S*.95),':','color','black','LineWidth',3);
hold on;
plot(ImpSort_pi(:,S*.5),'color','black','LineWidth',3);
hold on;
plot(ImpSort_pi(:,S*.05),':','color','black','LineWidth',3);
xlabel('Inflation response to productivity shock')

subplot(3,1,3);
hold on;
plot(ImpSort_y(:,S*.95),':','color','black','LineWidth',3);
hold on;
plot(ImpSort_y(:,S*.5),'color','black','LineWidth',3);
hold on;
plot(ImpSort_y(:,S*.05),':','color','black','LineWidth',3);
xlabel('Output response to productivity shock')



figure(2)
subplot(3,1,1);
plot(IMPR_r);
xlabel('Interest rate response to productivity shock')
subplot(3,1,2)
plot(IMPR_pi);
xlabel('Inflation response to productivity shock')
subplot(3,1,3)
plot(IMPR_y);
xlabel('Output response to productivity shock')

figure(3)
subplot(1,3,1);
hist(VARDECOMP(1,:),50);
subplot(1,3,2);
hist(VARDECOMP(2,:),50);
subplot(1,3,3);
hist(VARDECOMP(3,:),50);

figure(4)
plot_vardecomp(VARDECOMP,0)

figure(5)
for ss=1:9
    subplot(3,3,ss);
    hist(ImpSort_y(ss,:),50);
end


display('The posterior probability that output increases more than half a percent after a 1 s.d. productivity shocks is')
count/S

forec_sort=sort(forecast,3);
forec_low=reshape(forec_sort(:,:,S*0.05),3,periods+1);
forec_med=reshape(forec_sort(:,:,S*0.5),3,periods+1);
forec_upp=reshape(forec_sort(:,:,S*0.95),3,periods+1);

FC_low=[ NaN(3,T-1) forec_low;];
FC_med=[ NaN(3,T-1) forec_med;];
FC_upp=[ NaN(3,T-1) forec_upp;];
DATA=[Z NaN(3,periods);];

figure(6)
subplot(3,1,1);
plot(DATA(1,:))
hold on
plot(FC_low(1,:),':')
hold on
plot(FC_med(1,:),'--')
hold on
plot(FC_upp(1,:),':')
xlabel('Interest rate')

subplot(3,1,2);
plot(DATA(2,:))
hold on
plot(FC_low(2,:),':')
hold on
plot(FC_med(2,:),'--')
hold on
plot(FC_upp(2,:),':')
xlabel('Inflation')

subplot(3,1,3);
plot(DATA(3,:))
hold on
plot(FC_low(3,:),':')
hold on
plot(FC_med(3,:),'--')
hold on
plot(FC_upp(3,:),':')
xlabel('Output')

