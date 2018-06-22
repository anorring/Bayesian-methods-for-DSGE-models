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

%set paramters for artificial sample
T=100;%Sample length
S=25;% Cross-sectional survey response dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set starting values for parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example: man bites dog
rho1=0.9; %persistence of technology
rho2=0.7; %persistence of demand
sigu=.018;  %s.d. of state innov
sigud=0.01; %s.d. "demand" shock
sigur=.014;  %s.d. of m.p. shock
sigaj=.04;  %s.d. of island tech
sigzj1=.1;  %s.d. of private info noise
sigzj2=.1;  %s.d. of private info noise
sigdj=0.05;  %s.d. of island demand
sigmbd=0.1; %s.d. of m-b-d signal
varphi=1.1; %labour supply curvature
delta=1.07; %elasticity of demand
fir=0.08;%Interest inertia
fipi=1.55; %Taylor param;
fiy=0.24; %Taylor rule param
stick=0.74; %Calvo parameter
beta=0.99; %discount rate

theta=[rho1,rho2,sigu,sigud,sigur,sigaj,sigzj1,sigzj2,sigdj,sigmbd,varphi,delta,fir,fipi,fiy,stick,beta;]';
xmax=theta;
if startfrompreviousmax ==1;
    load('xmax');
    theta=xmax;
end
lastMCMC=[];

if restart==0;
    load('xmax');
    load('postmax');
    load('lldraw');
    load('priordraw');
    load('lastMCMC');
    theta=lastMCMC(:,end);
    load('vscale');
else
    postmax = -9e+200;
    vscale=diag(abs(theta))*0.00001+eye(length(theta));
end
%Define boundary of THETA
LB=[zeros(1,5),0,0,0,0,0,0,0,0,1.05,0,0,0.97;]';
UB=[0.9999,0.9999,ones(1,10)*5,1,3,2,1,1;]';
%-----------------------------------------------------------------
%Define hyper parameters etc
%-----------------------------------------------------------------
kbar=8;
tol=1e-4;


x=theta;


llcan = -9e+200;
priorcan=-9e+200;

%--------------------------------------------------------------------------
% Set control parameters
%--------------------------------------------------------------------------


%Number of draws (the ones in the outer loop of the MH)
S=1e4;
%Number of draws (the ones in the inner loop of the MH)
s=1e2;
%Calibrates vcscale for the first iterinitial iterations of S using an
%identity matrix times epseye.
iterinitial=13;
epseye=1e-5;
%Values for the constants in the adaptive proposal.
sd=5e-3; %Governs the acceptance rate : AIM FOR approx 0.23 !!
%Proportion of draws to be discarded for burnin purposes
burnin=0.5;
%Every how many iterations the program display the results
dispiter=100;
%Every how many iterations the program updates the scale matrix
calibiter=100;
%Every how many iterations the MCMC is saved
saveMCMC=1000;


fails=[];

%--------------------------------------------------------------------------
% Initializes the MH algorithm
%--------------------------------------------------------------------------

%Set up so that the first candidate draw is always accepted
if restart ==0
    lldraw = -9e+200;
    llcan = -9e+200;
    priordraw=-9e+200;
    priorcan=-9e+200;
    bdraw=x;
    mean0=x;
else
    
    bdraw=theta;
    [P,p,K,D,L,R,Rj,RRj,SigJ,M,N,a,b,dimx,dimX,dimu,dimuj,e1,e2,H,EE]= Lorenz(bdraw,kbar,tol); %#ok<NASGU>
    lpostdraw= MBDLL_Lorenz(M,N,a,b,e1,dimX,SigJ,bdraw,Z,T,CPImean,NGDPmean,NonZeroCPI,NonZeroNGDP,maxsurv,H);
    lpostcan = lpostdraw;
    
end



if restart==0;
    vscale=diag(abs(theta))*epseye;
    
end
%Store all draws in the following matrices which are initialized here
bb_=zeros(length(theta),s);
OutsideProp=zeros(S,1);
SwitchesProp=zeros(S,1);
%Number of draws outside parameter boundaries
q=0;
%Number of switches (acceptances)
pswitch=0;

%%
tic
[P,p,K,D,L,R,Rj,RRj,SigJ,M,N,a,b,dimx,dimX,dimu,dimuj,e1,e2,H,EE]= Lorenz(theta,kbar,tol);
Pdraw=P;pdraw=p;Mdraw=M;Kdraw=K;Ddraw=D;Ndraw=N;Ldraw=L;Rdraw=R;Rjdraw=Rj;RRjdraw=RRj;a_draw=a;b_draw=b;SigJdraw=SigJ;
EE %#ok<NOPTS>
toc
%%
%--------------------------------------------------------------------------
% MH algorithm starts here
%--------------------------------------------------------------------------

% Start of the outer MH loop.

for iter=1:S
    tic
    
    for iter2=1:s
        bcan = bdraw + norm_rnd(vscale);
        if min(bcan > LB)==1
            if min(bcan < UB)==1;
                [P,p,K,D,L,R,Rj,RRj,SigJ,M,N,a,b,dimx,dimX,dimu,dimuj,e1,e2,H,EE]= Lorenz(bcan,kbar,tol); 
                if EE==1;
                    priorcan=logpriorMBD(bcan);
                    llcan= MBDLL_Lorenz(M,N,a,b,e1,dimX,SigJ,bcan,Z,T,CPImean,NGDPmean,NonZeroCPI,NonZeroNGDP,maxsurv,H);
                    laccprob = llcan - lldraw + priorcan - priordraw;
                else
                    llcan= -9e+200;
                    laccprob= -9e+200;
                end
            else
                %                 counthis=[counthis 0];
                laccprob=-9e+200;
                q=q+1;
            end
        else
            laccprob=-9e+200;
            q=q+1;
        end
        if llcan + priorcan >= postmax;
            postmax=llcan + priordraw;
            xmax=bcan;
        end
        %Accept candidate draw with log prob = laccprob, else keep old draw
        if log(rand)<laccprob
            lldraw=llcan;
            priordraw=priorcan;
            bdraw=bcan;
            pswitch=pswitch+1;
            Pdraw=P;pdraw=p;Mdraw=M;Kdraw=K;Ddraw=D;Ndraw=N;Ldraw=L;Rdraw=R;Rjdraw=Rj;RRjdraw=RRj;a_draw=a;b_draw=b;SigJdraw=SigJ;
            
        end
        
        bb_(:,iter2)=bdraw;
        
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    toc
    
    disp('acceptance rate and fraction outside bounds');
    [pswitch/(iter*s) q/(iter*s)] %#ok<NOPTS>
    
    lastMCMC=[lastMCMC bb_(:,iter2);]; %#ok<AGROW>
    
    %     if length(lastMCMC(1,:)) >= (length(x)+1);
    if length(lastMCMC(1,:)) >= 100;
        SIGMCMC=cov(lastMCMC(:,50:end)');
        vscale=sd*SIGMCMC;
        if length(lastMCMC(1,:)) >= (length(x)+101);
            SIGold=cov(lastMCMC(:,1:end-100)');
            oldMCMC=lastMCMC(:,1:end-100);
        else
            SIGold=SIGMCMC;
            oldMCMC=lastMCMC;
        end
        disp('Current, mode, mean, change in mean last 10 000 draws, s.d., change in s.d. of MCMC last 10 000 draws');
        [bdraw, xmax, mean(lastMCMC')', (mean(lastMCMC')'-mean(oldMCMC')'), diag(SIGMCMC).^.5, (diag(SIGMCMC).^.5)-(diag(SIGold).^.5);] %#ok<NOPTS,UDIM>
        disp('Current posterior and at mode');[lldraw+priordraw postmax] %#ok<NOPTS>
        disp('Current likelihood and at mode');[lldraw postmax-logpriorMBD(xmax)] %#ok<NOPTS>
    else
        disp('Current mode');xmax %#ok<NOPTS>
        disp('Posterior Likelihood at mode');postmax %#ok<NOPTS>
        
    end
    disp('Total number of draws');length(lastMCMC(1,:))*s %#ok<NOPTS>
    
    
    %         end
    save('lastMCMC','lastMCMC');
    save('xmax','xmax');
    save('postmax','postmax');
    save('vscale','vscale');
    save('lldraw','lldraw');
    save('priordraw','priordraw');
    MCMC_Lorenz_25_prior=lastMCMC;
    save('MCMC_Lorenz_25_prior','MCMC_Lorenz_25_prior');
end

