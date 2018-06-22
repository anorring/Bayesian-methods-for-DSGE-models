%Program to plot probability intervals for IRFs etc using MCMC from M-H
%chain.
clc
clear all
close all
global Z
load('Z');
load('bb_');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=300;% NUmber of independent draws from the markov chain
periods=25; %Number of periods in the IRFs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find and plot probability intervals for IRFS, Variance decompisition and 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ra = max(size(bb_));
ff = ceil(ra.*rand(S,1));

IMPR=zeros(periods,S);%This is where we store the IRFs for the different draws
VARDECOMP=zeros(2,S);%This is where we store the variance decompositions for the different draws
for s=1:S;
    theta=bb_(:,ff(s,1));%Take a random draw from the MCMC
    [A,C,D,R]=DSGE_SS(theta);% Compute the state space matrices for theta
   %%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %Impulse response functiom
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
    IR=zeros(periods,1);
    for ss=0:periods-1
        IR(ss+1)= D(1,1)*(A^ss)*C(1,1);
    end
    IMPR(:,s)=IR;
    %%%%%%%%%%%%%%%%%%%%%%%%
    %Variance decomp
    %%%%%%%%%%%%%%%%%%%%%%%%
%     pivar=D(2)^2*inv(1-A^2)*sigx^2 + sigp^2;
%     VARDECOMP(:,s)=[(sigp^2)/pivar (D(2)^2)*(inv(1-r^2)*sigx^2)/pivar ]';
%     
    
end;

VARDECOMPsort=sort(VARDECOMP,2);
ImpSort=sort(IMPR,2);
figure
hold on;
plot(ImpSort(:,S*.95),':','color','black','LineWidth',3);
hold on;
plot(ImpSort(:,S*.5),'color','black','LineWidth',3);
hold on;
plot(ImpSort(:,S*.05),':','color','black','LineWidth',3);

figure
plot(IMPR)

figure
subplot(1,2,1);
hist(VARDECOMP(1,:),50);
subplot(1,2,2);
hist(VARDECOMP(2,:),50);
% close all
plotpost(VARDECOMP,0)