%DHOE_sim_and_estimate, simulates the data
clc
clear var;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign benchmark values to parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=100;%number of periods
S=25;%number of survey responses
kbar=40;%maximum order of expectations
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve model
[M,N,a,p_f_sd_j,H,Err]= DHOE_solve(theta,kbar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulate model
shocks=randn(3+S,T);
% surf(cov(shocks'))
X=zeros(kbar,T+1);
Z_p=zeros(1,T);
Z_f=zeros(S,T);
for t=2:T
    X(:,t)= M(1:kbar,1:kbar)*X(:,t-1)+N(1:kbar,:)*shocks(1:3,t);
    Z_p(1,t)=a(1:kbar)*X(:,t)+s_eps*shocks(2,t);     
    Z_f(1:S,t)=ones(S,1)*a(1:kbar)*M(1:kbar,1:kbar)*H(1:kbar,1:kbar)*X(:,t) + p_f_sd_j*shocks(4:end,t);
end
figure
subplot(2,1,1);
plot(Z_p,'linewidth',2)
xlabel('Price of asset')
subplot(2,1,2);
plot(Z_f')
xlabel('One period ahead price forecasts')

save Z_p Z_p
save Z_f Z_f
[LL]= DHOE_LL(theta,kbar,Z_p,Z_f)
