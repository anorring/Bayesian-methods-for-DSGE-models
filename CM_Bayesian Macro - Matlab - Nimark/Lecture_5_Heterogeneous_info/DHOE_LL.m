function [LL]= DHOE_LL(theta,kbar,Z_p,Z_f)

[M,N,a,p_f_sd_j,H,Err]= DHOE_solve(theta,kbar);
if Err==0
T=length(Z_p);
S=length(Z_f(:,1));
A=M;
C=[N zeros(kbar+1,S);];

D=[a;
    ones(S,1)*a(1:end-1)*M(1:end-1,1:end-1)*H;];

R=[0,1,0,zeros(1,S);
    zeros(S,3), eye(S)*p_f_sd_j;];

Z=[Z_p;
    Z_f;];
Xfilt=zeros*M(:,1);
P=dlyap(M,N*N');
LL=0;
for t=1:T
    Z_tilda=Z(:,t)-D*Xfilt;%These are the innovations (i.e. Ztilde)
    Omega=D*A*P*A'*D'+(D*C+R)*(D*C+R)';
    Omegainv=eye(S+1)/Omega;
    K=(A*P*(D*A)'+C*C'*D'+C*R')*Omegainv;
    P=A*P*A'+C*C'-K*Omega*K';
    Xfilt=A*Xfilt+A*K*Z_tilda;
    P = A*(P-P*D'*Omegainv*D*P)*A' + C*C';
    LL = LL - 0.5*(log(det(Omega)) + Z_tilda'*Omegainv*Z_tilda);
end
if imag(LL)~=0
    LL=-9e+200;
end
else
    LL=-9e+200;
end