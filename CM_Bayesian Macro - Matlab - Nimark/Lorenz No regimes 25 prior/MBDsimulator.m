% Man-Bites-Dog simulator

clc
clear all
close all

rho=0.95; %persistence of state
sigu=1;%S.d. of state innov
sigeta=1;%s.d. of private info noise
sigeps=1;%s.d. of public signal noise
alfa=0.1;%unconditional prob of Z=1, i.e. of observing pub signal
gamma=3;%s.d. multiplier of u when Z=1 

tol=1e-4; %convergence criteria
bindim=6;% Number of periods that matter (Markov order, if you like) and dimension of binary identifier of time varying matrices
binmax=2^(bindim); % Number of "regimes"
for j=1:bindim;
    binbase(1,bindim-j+1)=2^(j-1);
end

theta=[rho,sigu,sigeta,sigeps,alfa,gamma]';
% load('xmax');
% theta=xmax;
%MBDplotter(theta);
%%
% [M,N,A,DD0,DD1,KK0,KK1,R0,R1,P,p]=MBDsolve(theta);
[M,N,A,DD0,DD1,KK0,KK1,R0,R1,P,p]=MBDsolve(theta);

T=100+bindim;Z=zeros(1,T);
for t=1:T;
    if (rand) <= alfa ;
        Z(1,t)=1;
    end
end

U=randn(3,T);
X=zeros(size(M,1),T);
S=zeros(1,T);

for t=bindim:T
    z=Z(1,t-bindim+1:t);
    j=binbase*z'+1;
    
    X(:,t)= M(:,:,j)*X(:,t-1)+N(:,:,j)*U(:,t);
    S(t)=A(:,:,j)*X(:,t);
    
end


Xc=zeros(size(M,1),T);
Sc=zeros(1,T);
Uc=U;
for t=1:T;
    if Z(t)==1;
        Uc(1,t)=theta(6)*U(1,t);
    end
end

Z0=Z*0;
for t=bindim:T
    z=Z0(1,t-bindim+1:t);
    j=binbase*z'+1;    
    Xc(:,t)= M(:,:,j)*X(:,t-1)+N(:,:,j)*Uc(:,t);
    Sc(t)=A(:,:,j)*Xc(:,t);
    
end
    
   figure
   subplot(2,1,1);
   plot(S(5:end),'r'); hold on; plot(Sc(5:end));
   subplot(2,1,2);
   plot(Z(5:end));

        
        

