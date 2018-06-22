function [P,p,K,D,L,R,Rj,RRj,SigJ,M,N,a,b,dimx,dimX,dimu,dimuj,e1,e2,H,EE]= Lorenz(theta,kbar,tol)
EE=1;
%%
try
    
rho1=theta(1); %persistence of technology
rho2=theta(2); %persistence of demand
sigu=theta(3);  %s.d. of state innov
sigud=theta(4); %s.d. "demand" shock
sigur=theta(5);  %s.d. of m.p. shock
sigaj=theta(6);  %s.d. of island tech
sigzj1=theta(7);  %s.d. of private info noise
sigzj2=theta(8);  %s.d. of private info noise
sigdj=theta(9);  %s.d. of private info noise
sigmbd=theta(10); %s.d. of m-b-d signal
varphi=theta(11); %labour supply curvature
delta=theta(12); %elasticity of demand
fir=theta(13);%Interest inertia
fipi=(1-fir)*theta(14); %Taylor param;
fiy=(1-fir)*theta(15); %Taylor rule param
stick=theta(16); %Calvo parameter
beta=theta(17); %discount rate

lambda=(1-stick)*(1-stick*beta)/stick;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define dimensions etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimx=2;
dimX=dimx*kbar;
dimZ=6;
dimuj=4;
dimu=4;


if dimZ >= dimu+dimuj;
    disp('Dear Sir, you may have a perfectly revealing equilibrium.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define matrix arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exogenously specified matrices
C=zeros(dimx,dimu);
R=zeros(dimZ,dimu);
Rj=zeros(dimZ,dimuj);

%Other useful stuff
e1=[1,zeros(1,dimX-1)];
e2=[0,1,zeros(1,dimX-2)];
H=[zeros(dimX,dimx) eye(dimX);];
H=H(:,1:end-dimx);

%Innovations to aggregate productivity and demand
C(1,1)=sigu;
C(2,2)=sigud;

%Private noise
Rj(1,1)=sigaj;
Rj(1,1)=sigaj;
Rj(2,4)=sigdj;
Rj(3,2)=sigzj1;
Rj(4,3)=sigzj2;

%Public noise
R(5,:)=[0,0,0,sigur;];
R(6,:)=[0,0,sigmbd,0;];

D=zeros(dimZ,dimX); %D used to find starting values
D(1:2,1:2)=eye(2);
D(3:6,1)=1;

% P=zeros(dimX,dimX);
p=zeros(dimX,dimX);
K=zeros(dimX,dimZ);
% L=zeros(dimZ,dimZ);
a=zeros(1,dimX);
b=zeros(1,dimX);
EE=1;




%Endogenously specified matrices
% D=zeros(dimZ,dimX);
Mst=zeros(dimX,dimX);
Nst=zeros(dimX,dimu+dimuj);

Ast=[rho1,0;0,rho2;];
for j=1:dimx:dimX;
    Mst(j:j+dimx-1,j:j+dimx-1)=Ast;
    Nst(j:j+dimx-1,1:dimu)=C;
end
N=Nst;

ast=(lambda -lambda*delta*varphi)*e1;
bst=(ast*Mst*H)-fipi*ast;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define and assign hyperparamters for iterative solution loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


maxiter=2000;
M=Mst;
Nst=N;

RRj=[R Rj];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Step=ones(4,1)*0.1;
Hstep=0;
Diff=ones(1,4);DiffOld=ones(1,4);iter=0;


p=eye(dimX);

P=eye(dimX);

SigJ=eye(dimX);
L=eye(dimZ);

%%
while max(abs(Diff)) > tol && iter < maxiter;
    %%
  if max(max(abs(M*H)))<= 1;
    a=lambda*(bst-e1)+lambda*delta*varphi*(bst-e1) + beta*ast*Mst*H;
    b=e2+(ast+bst)*Mst*H-fipi*ast-fiy*bst;
  else
      a=ast;
      b=bst;
  end 

        D(3,:)=a;
        D(4,:)=delta*a + b;
        D(5,:)=(1-fir)*fipi*a + (1-fir)*fiy*b;   
        
 
    %The Nimark Filter option
    L=(D*M)*p*(D*M)'+(D*N+RRj)*(D*N+RRj)';
    K=(M*p*(D*M)'+N*N'*D'+N*RRj')/(L);
    p=P-K*L*K';
    P=M*p*M'+N*N';
    
    
    KDM=K*D*M;
    M(dimx+1:end,:)=[KDM(1:end-dimx,1:end-dimx) zeros(dimX-dimx,dimx)] + [zeros(dimX-dimx,dimx) M(1:end-dimx,1:end-dimx)] - [zeros(dimX-dimx,dimx) KDM(1:end-dimx,1:end-dimx)];
    
    
    KDN=K*D*N;
    KR=K*[R Rj*0];
    N(dimx+1:dimX,:)=KDN(1:end-dimx,:) + KR(1:end-dimx,:) ;
    
   
    SigJ=(eye(dimX)-K*D)*M*SigJ*M'*(eye(dimX)-K*D)'+K*Rj*Rj'*K';
    
    Step=(1-Hstep)*Step+Hstep*Step.*(DiffOld./Diff)';
    Step(1)=max([Step(1),0.00001;]);
    Step(2)=max([Step(2),0.00001;]);
    Step(3)=max([Step(3),0.1;]);
    Step(4)=max([Step(4),0.1;]);
    
    
    
    Step(1)=min([Step(1),1;]);
    Step(2)=min([Step(2),1;]);
    Step(3)=min([Step(3),1;]);
    Step(4)=min([Step(4),1;]);
    
    
    DiffOld=Diff;
    DiffM=max(max(abs(M-Mst)));
    DiffN=max(max(abs(N-Nst)));
    Diffa=max(max(abs(a-ast)));
    Diffb=max(max(abs(b-bst)));
    
    Diff=[DiffM,DiffN,Diffa,Diffb;];
    
    Mst=Step(1)*M+(1-Step(1))*Mst;M=Mst;
    Nst=Step(2)*N+(1-Step(2))*Nst;N=Nst;
    ast=Step(3)*a+(1-Step(3))*ast;a=ast;
    bst=Step(4)*b+(1-Step(4))*bst;b=bst;
   
    
    iter=iter+1;
  
    if iter > 2000;
        EE=0;
        
        M=[];
        break
    end
end

catch
    EE=0;
    M=[];
end