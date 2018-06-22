%Finds maximum likelihood estiamte of parameters for univariate UC model
%with unit root
% 	tau[t] = tau[t-1] + Cu[t]
% 
%   pi[t] = tau[t] + Ru[t]
%
clear all;
clc
close all

load USdata

Z=USdata(1,:);

%Scalar UC values for A and C in state space system
A=1;
D=1;

dimX=size(A,1);
dimZ=length(D(:,1));

% periods=1000;
randn('seed', 12345);

Cgrid=[0.0001:0.0001:0.01];
Rgrid=[0.0001:0.001:0.01];
dimC=length(Cgrid);
dimR=length(Rgrid);

T=length(Z);

Grid=zeros(dimC,dimR);

for c=1:dimC;
 
    for r=1:dimR;
        C=Cgrid(c);
        R=Rgrid(r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Kalman filter equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=zeros(dimX,dimX,T+1);
K=zeros(dimX,dimZ,T);
Xtt=zeros(1,T+1);
LL=0;
P(:,:,1)=C*C'/(1-0.9*A^2);

for t=1:T;
    P(:,:,t+1)=A*(P(:,:,t)-(P(:,:,t)*D'/(D*P(:,:,t)*D'+R*R'))*D*P(:,:,t))*A'+C*C';    
    K(:,:,t)= P(:,:,t+1)*D'/(D*P(:,:,t+1)*D'+R*R');
    Xtt(:,t+1)=A*Xtt(:,t)+K(:,:,t)*(Z(:,t)-D*A*Xtt(:,t));
    Ztilda=(Z(:,t)-D*A*Xtt(:,t));
    Omega=D*P(:,:,t+1)*D'+R*R';
    LL=LL-0.5*(log(2*pi)+log(det(Omega))+(Ztilda'/Omega)*Ztilda);
end
Grid(c,r)=LL;
    end
end
[LLmaxRow IRow]=max(Grid);
[LLmax I]=max(LLmaxRow);
maxindex=[IRow(I) I];
MLEtheta=[Cgrid(maxindex(1)) Rgrid(maxindex(2))]
C=Cgrid(maxindex(1));
R=Rgrid(maxindex(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Kalman filter equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=zeros(dimX,dimX,T+1);
K=zeros(dimX,dimZ,T);
Xtt=zeros(1,T+1);
P(:,:,1)=C*C'/(1-0.9*A^2);
for t=1:T;
    P(:,:,t+1)=A*(P(:,:,t)-(P(:,:,t)*D'/(D*P(:,:,t)*D'+R*R'))*D*P(:,:,t))*A'+C*C';    
    K(:,:,t)= P(:,:,t+1)*D'/(D*P(:,:,t+1)*D'+R*R');
    Xtt(:,t+1)=A*Xtt(:,t)+K(:,:,t)*(Z(:,t)-D*A*Xtt(:,t));
end

figure(1)

plot(Z,'linewidth',2)
hold on;
plot(Xtt(2:end),'linewidth',2,'linestyle','--','color','r')
legend('\pi','\tau_{t|t}')
CC=[C,0;];
RR=[0,R;];
[Xsmooth]=smooth(A,CC,D,RR,Z);
figure(2)
plot(Xtt(2:end),'linewidth',2,'linestyle','-')
hold on
plot(Xsmooth(2:end),'linewidth',2,'linestyle','--','color','g')
legend('\tau_{t|t}','\tau_{t|T}')

[M,L,U]=distplot(A,CC,D,RR,Z,0.95,0.05,100,1,0);

