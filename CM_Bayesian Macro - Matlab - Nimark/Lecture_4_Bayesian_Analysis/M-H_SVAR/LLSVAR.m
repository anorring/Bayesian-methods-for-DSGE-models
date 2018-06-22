function [LL]=LLSVAR(theta,Z)

C=[theta(1),0;
   theta(2),theta(3);];

A=[theta(4),theta(5);
    theta(6),theta(7);];

CC=C*C';

RES=Z(:,2:end)-A*Z(:,1:end-1);
T=length(Z); 

    
 LL=-(T*2/2)*log(2*pi)-(T/2)*log(det(CC))-0.5*trace((RES'/CC)*RES);
if max(abs(eig(A))) >= 1;
    ll=-1d20;
end