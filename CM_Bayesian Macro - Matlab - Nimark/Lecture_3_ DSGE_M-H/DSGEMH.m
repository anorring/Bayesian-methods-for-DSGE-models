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
J=1e4;
epseye=1e-4;

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
            lpostcan = log_prior_DSGE(bcan)+LLDSGE(bcan);%switch on for use of priors
%               lpostcan = log_prior_DSGE(bcan);%Uncomment for prior predictive analysis
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
        vscale=5e-4*cov(bb_(:,1:iter)');
        iter
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

convcheck(bb_(:,100:end));

% figure
plotpost(bb_,0)
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Find and plot probability intervals for IRFS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% periods=25;
% ra = max(size(bb_));
% 
% S=300;
% ff = ceil(ra.*rand(S,1));
% 
% IMPR=zeros(periods,S);
% VARDECOMP=zeros(2,S);
% for s=1:S;
%     theta=bb_(:,ff(s,1));
%     r=theta(1); %productivity persistence
%     g=theta(2); %relative risk aversion
%     d=theta(3); %Calvo parameter
%     b=theta(4);   %discount factor
%     k=((1-d)*(1-d*b))/d; %slope of Phillips curve
%     f=theta(5);% coefficient on inflation in Taylor rule
%     sigx=theta(6);% s.d. prod shock
%     sigy=theta(7);% s.d. demand shock
%     sigp=theta(8);% s.d. cost push shock
%     c=g-k*r-2*g*r+k*f+g*(r^2);
%     SIGvv=diag([sigy;sigp;]);
%     A=r;
%     C=sigx;
%     D=[(-k*(f-r))/-c;(k*(1-r))/-c;];
%     
%     IR=zeros(periods,1);
%     for ss=0:periods-1
%         IR(ss+1)= f*D(2,1)*(r^ss)*C;
%     end
%     IMPR(:,s)=IR;
%     %Variance decomp
%     pivar=D(2)^2*inv(1-r^2)*sigx^2 + sigp^2;
%     VARDECOMP(:,s)=[(sigp^2)/pivar (D(2)^2)*(inv(1-r^2)*sigx^2)/pivar ]';
%     
% end;
% VARDECOMPsort=sort(VARDECOMP,2);
% ImpSort=sort(IMPR,2);
% figure
% hold on;
% plot(ImpSort(:,S*.95),':','color','black','LineWidth',3);
% hold on;
% plot(ImpSort(:,S*.5),'color','black','LineWidth',3);
% hold on;
% plot(ImpSort(:,S*.05),':','color','black','LineWidth',3);
% 
% figure
% plot(IMPR)
% 
% figure
% subplot(1,2,1);
% hist(VARDECOMP(1,:),50);
% subplot(1,2,2);
% hist(VARDECOMP(2,:),50);
% % close all
% plotpost(VARDECOMP,0)
