function [LL]=LLDSGE(theta)
global Z

r=theta(1); %productivity persistence
g=theta(2); %relative risk aversion
d=theta(3); %Calvo parameter
b=theta(4);   %discount factor
k=((1-d)*(1-d*b))/d; %slope of Phillips curve
f=theta(5);% coefficient on inflation in Taylor rule
sigx=theta(6);% s.d. prod shock
sigy=theta(7);% s.d. demand shock
sigp=theta(8);% s.d. cost push shock
sigr=theta(9);% s.d. cost push shock

c=g-k*r-2*g*r+k*f+g*(r^2); 
  
A1=[1,-f,0;
    0,1,-k;
    (1/g),0,1;];
R1=[sigr,0,0;0,sigp,0;0,0,sigy;];
R=A1\R1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put model in state space form
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=r;
C=sigx;
D=[f*((k*(1-r))/-c);(k*(1-r))/-c;(-k*g*(f-r))/-c;];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial values for Kalman filter etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n=length(D(:,1)); %dimension of observables;    
    T=length(Z);       
    LL=0;    
    P= inv(1-r^2)*sigx^2;%Initial uncertainty equal to unconditional variance of state
    Xfilt=0; %initial value for the filtered state.
%    RR=SIGvv*SIGvv';
RR=R*R';
   QQ=C*C';
   
       
    %Compute recursive likleihood using the Kalman filter 
    for tt=1:T
        a=Z(:,tt)-D*Xfilt;%These are the innovations (i.e. Ztilde)
        Omega=D*P*D'+RR;
        Omegainv=inv(Omega);
        K=P*D'*Omegainv;
        Xfilt=A*Xfilt+A*K*a;
        P = A*(P-P*D'*Omegainv*D*P)*A' + QQ;
        LL = LL - 0.5*(log(det(Omega)) + a'*Omegainv*a);
    end
    if imag(LL)~=0
        LL=-9e+200;
    end
    
    LL=-LL;%Because we are minimizing
