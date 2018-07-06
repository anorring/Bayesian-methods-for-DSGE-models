%solves model in DHOE
function [M,N,a,p_f_sd_j,H,Err]= DHOE_solve(theta,kbar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign benchmark values to parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=theta(1);%discount factor
r=theta(2);%persistence of theta
s_u=theta(3);%s.d. persistent shock
s_eps=theta(4);%s.d. transitory shock
s_eta=theta(5);%idisoyncratic measurement error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define some useful matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R1=[0,0;0,-s_eps;];
R2=[s_eta;0;];
%R=[R1 R2];
%dimj=1;%dimension of private signal vector

Err=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values for k=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=r;
N=[s_u,0,0;];
a=-1;
A=NaN(kbar,kbar);
A(1,1)=a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main solution loop with increasing dimension k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    for k=1:kbar
        
        e1=[1,zeros(1,k-1);];
        L=[e1;a;];
        [K,~]=Kalman(M,N,L,R1,R2);
        
        Mst=[r,zeros(1,k);zeros(k,1),zeros(k,k);] + [zeros(1,k),zeros(1,1);K*L*M,zeros(k,1);] + [zeros(1,1),zeros(1,k);zeros(k,1),(eye(k)-K*L)*M;];
        Nst=[s_u,0,0;
            K*L*N+K*[R1 0*R2];];
        H=[zeros(k,1),eye(k);];
        a=[-e1,0;]+b*a*M*H;
        A(k+1,1:k+1)=a;
        M=Mst;
        N=Nst;
    end
    IKLM=(eye(kbar)-K*L)*M(1:kbar,1:kbar);
    KRRK=K*(R2*R2')*K';
    SIGj=dlyap(IKLM,KRRK);
    p_f_sd_j=((a(1:end-1)*M(1:end-1,1:end-1))*SIGj*(a(1:end-1)*M(1:end-1,1:end-1))')^0.5;
    %
    if max(abs(eig(M)))>=1
        Err=1;
    end
catch
    Err=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

