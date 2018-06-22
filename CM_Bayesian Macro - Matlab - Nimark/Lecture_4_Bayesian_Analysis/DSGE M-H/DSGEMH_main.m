% Set up and estimate miniture DSGE model
clc
clear all
close all
global Z
load('Z');
%--------------------------------------------------------------------------
% Set control parameters
%--------------------------------------------------------------------------

%Number of draws
J=5e5;
epseye=1e-3;
burnin=J*0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define structural parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=0.95; %productivity persistence
g=3; %relative risk aversion
d=0.75; %Calvo parameter
b=0.99;   %discount factor
k=((1-d)*(1-d*b))/d; %slope of Phillips curve
f=1.5;% coefficient on inflation in Taylor rule
sigx=0.1;% s.d. prod shock
sigy=0.1;% s.d. demand shock
sigp=0.1;% s.d. cost push shock
sigr=0.1;% s.d. cost push shock


theta=[r,g,d,b,f,sigx,sigy,sigp,sigr]';%Starting value for parameter vector
LB=[0,1,0,0,1,zeros(1,4);]';
UB=[1,10,1,1,5,1*ones(1,4);]';
x=theta;



%--------------------------------------------------------------------------
% Initializes the MH algorithm
%--------------------------------------------------------------------------

%Initializes the proposal variance.
%Set up so that the first candidate draw is always accepted
lpostdraw = -9e+200;
bdraw=x;
vscale=diag(abs(theta))*epseye+1e-5*eye(length(x));


bb_=zeros(length(x),J);%Store all draws in bb_
%Matrices that keep track of switches and drwas outside LB and UB
OutsideProp=zeros(J,1);
SwitchesProp=zeros(J,1);
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
for iter=1:J
    iter=iter+1;
    % Draw from proposal density Theta*_{t+1} ~ N(Theta_{t},vscale)
    bcan = bdraw + norm_rnd(vscale);
    
    if min(bcan > LB)==1
        if min(bcan < UB)==1
%             lpostcan = LLDSGE(bcan);%Uncomment for improper uniform priors
%             lpostcan = log_prior_DSGE(bcan)+LLDSGE(bcan);%switch on for use of priors
              lpostcan = log_prior_DSGE(bcan);%Uncomment for prior predictive analysis
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
    
    bb_(:,iter)=bdraw;
    
    OutsideProp(iter)=q/iter;
    SwitchesProp(iter)=pswitch/iter;
    
    if iter >= 50 && mod(iter,1000)==0
        vscale=5e-1*cov(bb_(:,1:iter)');
        iter
        SwitchesProp(iter)
    end
    
end
toc

disp(['iter: ',num2str(iter)]);
disp(['acceptance rate: ',num2str(SwitchesProp(iter))]);

figure
bb_=bb_(:,50:end);
for j=1:8;
    subplot(3,3,j);
    hist(bb_(j,:),50);
end

convcheck(bb_(:,10:end));

% figure
plotpost(bb_(:,burnin:end),0)
load bb_mle
plotpost(bb_mle(:,burnin:end),1)
load bb_prior
plotpost(bb_prior(:,burnin:end),0)
%%

