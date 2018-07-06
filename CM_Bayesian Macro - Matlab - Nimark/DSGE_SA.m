% Set up and estimate miniture DSGE model
clc
clear all
close all
global Z
load('Z'); %load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial values of structural parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=0.95; %productivity persistence
g=5; %relative risk aversion (gamma)
d=0.75; %Calvo parameter (price stickyness)
b=0.99;   %discount factor
k=((1-d)*(1-d*b))/d; %slope of Phillips curve
f=1.5;% coefficient on inflation in Taylor rule (phi)
sigx=0.1;% s.d. prod shock
sigy=0.11;% s.d. demand shock
sigp=0.1;% s.d. cost push shock
sigr=0.1;% s.d. cost push shock

names={'r' 'g' 'd' 'b' 'f' 'sigx' 'sigy' 'sigp' 'sigr'};
theta=[r,g,d,b,f,sigx,sigy,sigp,sigr]';%Starting value for parameter vector
LB=[0,0,0,0,1,zeros(1,4);]';%Lower bound for parameter vector
UB=[1,10,1,1,5,1*ones(1,4);]';%Upper bound for parameter vector
x=theta;

% x=THETA;
% initial temperature
sa_t= 5;
% temperature reduction factor
sa_rt=.3;
% number of draws before the change in temperature
sa_nt=5;
% number of draws before a step size adjustment
sa_ns=5;
% warning off all;

[xhat]=simannb( 'LLDSGE', x, LB, UB, sa_t, sa_rt, sa_nt, sa_ns, 1);

%--------------------------------------------------------------------------

theta=xhat;
r=theta(1); %productivity persistence
g=theta(2); %relative risk aversion
d=theta(3); %Calvo parameter
b=theta(4);   %discount factor
k=((1-d)*(1-d*b))/d; %slope of Phillips curve
f=theta(5);% coefficient on inflation in Taylor rule
sigx=theta(6);% s.d. prod shock
sigy=theta(7);% s.d. demand shock
sigp=theta(8);% s.d. cost push shock
sigr=theta(9);% s.d. measurement error interest rate

c=g-k*r-2*g*r+k*f+g*(r^2); 
A1=[1,-f,0;
    0,1,-k;
    (1/g),0,1;];
R1=[sigr,0,0;0,sigp,0;0,0,sigy;];
SIGvv=[zeros(3,1), A1\R1;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put model in state space form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=r;
C=[sigx,zeros(1,3);];
D=[f*(k*(1-r))/-c;(k*(1-r))/-c;((-k*g*(f-r))/-c)];

[X]=smooth(A,C,D,0,SIGvv,Z);
sampleDates=[1981.5:.25:2010.5];
subplot(2,1,1,'fontsize',18);plot(sampleDates,Z','linewidth',2);
legend('FFR','Inflation','Output')
hold on;
subplot(2,1,2,'fontsize',18); plot(sampleDates,X(:,2:end),'linewidth',2);
legend('Smoothed Estimate of Productivity')

