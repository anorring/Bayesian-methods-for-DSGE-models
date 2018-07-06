%M-H for model from Dynamic Higher Order Expectations using simulated data
clc
clear var;
close all;
load Z_p
load Z_f
%plot data
figure
subplot(2,1,1);
plot(Z_p,'linewidth',2)
title('Price of asset')
subplot(2,1,2);
plot(Z_f')
title('One period ahead price forecasts')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set parameters for Assign benchmark values to parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=length(Z_f(:,1));%number of survey responses
kbar=40; %maximum order of expectations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign benchmark values to parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=0.9;%discount factor
r=0.9;%persistence of theta
s_u=.1;%s.d. persistent shock 
s_eps=.5;%s.d. transitory shock
s_eta=.2;%idisoyncratic measurement error
theta=[b,r,s_u,s_eps,s_eta]';

%set bounds of parameter space
LB=[0,0,0,0,0;];UB=[0.99,0.99,5,5,5;];


%--------------------------------------------------------------------------
% Initializes the MH algorithm
%--------------------------------------------------------------------------
J=1000;%Number of draws 
%Initializes the proposal variance.
%Set up so that the first candidate draw is always accepted
lpostdraw = -9e+200;
theta_draw=theta;
vscale=diag(1e-5*abs(theta)) + 1e-5*eye(length(theta));

MCMC=zeros(length(theta),J);%Store all draws in bb_

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
    theta_can = theta_draw + norm_rnd(vscale);    
    if min(theta_can > LB)==1
        if min(theta_can < UB)==1
            lpostcan = DHOE_LL(theta_can,kbar,Z_p,Z_f);
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
        theta_draw=theta_can;
        pswitch=pswitch+1;
    end
    
    MCMC(:,iter)=theta_draw;
    
    OutsideProp(iter)=q/iter;
    SwitchesProp(iter)=pswitch/iter;
    
    if iter >= 100 && mod(iter,100)==0
        vscale=1e-3*cov(MCMC(:,1:iter)')+ 1e-6*eye(length(theta));%Adaptive proposal density
        
disp(['iter: ',num2str(iter)]);
disp(['acceptance rate: ',num2str(SwitchesProp(iter))]);
        theta_draw
        lpostdraw
                
    end
    
end
toc


%plot raw MCMC
figure
for jj=1:5;
    subplot(2,3,jj);
    plot(MCMC(jj,:))
    if jj==1
        xlabel('\beta')
    end
    if jj==2
        xlabel('\rho')
    end
    if jj==3
        xlabel('\sigma_u')
    end
    if jj==4
        xlabel('\sigma_{\epsilon}')
    end
    if jj==5
        xlabel('\sigma_{\eta}')
    end
end

%plot marginal posterior parameter densities
figure
for jj=1:5;
    subplot(2,3,jj);
    
    [f,xi] = ksdensity(MCMC(jj,:)');
    plot(xi,f,'LineWidth',2)
    hold on
    if jj==1
        xlabel('\beta')
    end
    if jj==2
        xlabel('\rho')
    end
    if jj==3
        xlabel('\sigma_u')
    end
    if jj==4
        xlabel('\sigma_{\epsilon}')
    end
    if jj==5
        xlabel('\sigma_{\eta}')
    end
end






