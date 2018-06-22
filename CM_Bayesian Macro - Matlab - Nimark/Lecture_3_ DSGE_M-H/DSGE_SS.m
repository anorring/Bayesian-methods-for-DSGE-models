function [A,C,D,R]=DSGE_SS(theta)



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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put model in state space form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=r;
C=[sigx,0,0,0;];
D=[f*((k*(1-r))/-c);(k*(1-r))/-c;(-k*g*(f-r))/-c;];
A1=[1,-f,0;
    0,1,-k;
    (1/g),0,1;];
R1=[sigr,0,0;0,sigp,0;0,0,sigy;];
R=A1\R1;

R=[zeros(3,1),R;];
