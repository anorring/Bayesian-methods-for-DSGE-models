%--------------------------------------------------------------------------
% Random Walk Chain Metropolis-Hastings algorithm to estimate stationary
% version of Lorenzoni (2009) on simulated data.
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
format short;
warning off all;
kbar=10;
tol=1e-6;
%set parameters for artificial sample
T=100;%Sample length
S=25;% Cross-sectional survey response dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set starting values for parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example: man bites dog
rho1=0.9; %persistence of technology
rho2=0.7; %persistence of demand
sigu=.02;  %s.d. of state innov
sigud=0.01; %s.d. "demand" shock
sigur=.01;  %s.d. of m.p. shock
sigaj=.04;  %s.d. of island tech
sigzj1=.1;  %s.d. of private info noise
sigzj2=.1;  %s.d. of private info noise
sigdj=0.05;  %s.d. of island demand
sigps=0.01;  %s.d. of public signal noise
varphi=1.5; %labour supply curvature
delta=2; %elasticity of demand
fir=0.08;%Interest inertia
fipi=1.55; %Taylor param;
fiy=0.24; %Taylor rule param
stick=0.74; %Calvo parameter
beta=0.99; %discount rate

theta=[rho1,rho2,sigu,sigud,sigur,sigaj,sigzj1,sigzj2,sigdj,sigps,varphi,delta,fir,fipi,fiy,stick,beta;]';
%solve model
[P,p,K,D,L,R,Rj,RRj,SigJ,M,N,a,b,dimx,dimX,dimu,dimuj,e1,e2,H,EE]= Lorenz(theta,kbar,tol);

% simulate data
for t=1:T;
    
