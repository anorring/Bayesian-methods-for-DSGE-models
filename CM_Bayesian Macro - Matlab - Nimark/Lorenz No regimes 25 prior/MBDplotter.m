function MBDplotter(theta)


%Man Bites Dog simple example: Only exogenous signals (c) K Nimark 2010


[M,N,A,DD0,DD1,KK0,KK1,R0,R1,P,p]=MBDsolve(theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define hyper parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kbar=10;% Number of orders of expectations
tol=1e-4; %convergence criteria
bindim=6;% Number of periods that matter (Markov order, if you like) and dimension of binary identifier of time varying matrices
binmax=2^(bindim); % Number of "regimes"
for j=1:bindim;
    binbase(1,bindim-j+1)=2^(j-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=theta(1); %persistence of state
sigu=theta(2);%S.d. of state innov
sigeta=theta(3);%s.d. of private info noise
sigeps=theta(4);%s.d. of public signal noise
alfa=theta(5);%unconditional prob of Z=1, i.e. of observing pub signal
gamma=theta(6);%s.d. multiplier of u when Z=1 



H=[zeros(kbar,1), eye(kbar);];H=H(:,1:kbar);
e1=[1,zeros(1,kbar-1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price and hierarchy IRF
periods=20;
ZVEC=[zeros(1, bindim-1) 1 zeros(1,periods-1) ];
j=binbase*ZVEC(1,(1:bindim))';
if ZVEC(1,bindim)==0;
    state(:,1)=N(:,1,j+1);
else
state(:,1)=N(:,1,j+1);
end
for t=1:periods;
    j=binbase*ZVEC(1,t:t+bindim-1)';
    state(:,t+1)=M(:,:,j+1)*state(:,t);
    price(1,t)=A(1,:,j+1)*state(:,t);
end
figure(1)
subplot(3,1,1);plot(price);hold on;
subplot(3,1,2);plot(state');hold on;

ZVEC=[zeros(1, bindim-1) 0 zeros(1,periods-1) ];
j=binbase*ZVEC(1,(1:bindim))';
if ZVEC(1,bindim)==0;
    state(:,1)=N(:,1,j+1)*gamma;
else
state(:,1)=N(:,1,j+1)*gamma;
end
for t=1:periods;
    j=binbase*ZVEC(1,t:t+bindim-1)';
    state(:,t+1)=M(:,:,j+1)*state(:,t);
    price(1,t)=A(1,:,j+1)*state(:,t);
end
figure(1)
subplot(3,1,1);plot(price,'--');hold on;
subplot(3,1,3);plot(state');hold on;

% %IRF of cross-sectional dispersion

periods=20;
ZVEC=[zeros(1, bindim) 1 zeros(1,periods-1) ]; %define sequence of Z
j=binbase*ZVEC(1,(1:bindim))';

LM=(eye(kbar)-KK0(:,:,j+1)*DD0(:,:,j+1))*M(:,:,j+1);
ISV=dlyap(LM,KK0(:,:,j+1)*R0(1,:)*R0(1,:)'*KK0(:,:,j+1)');
DISP1=[];
DISP1(1,1)= e1*ISV*(e1)';
SIG=ISV;
for t=1:periods;
    
    j=binbase*ZVEC(1,t+1:t+bindim)';
    
    if ZVEC(1,t+bindim)==0;
        LM=(eye(kbar)-KK0(:,:,j+1)*DD0(:,:,j+1))*M(:,:,j+1);
        KRRK0=KK0(:,:,j+1)*R0*R0'*KK0(:,:,j+1)';
        SIG=LM*SIG*LM'+ KRRK0;

     
    else
       LM=(eye(kbar)-KK1(:,:,j+1)*DD1(:,:,j+1))*M(:,:,j+1);
       KRRK1=KK1(:,:,j+1)*R1*R1'*KK1(:,:,j+1)';
       SIG=LM*SIG*LM'+ KRRK1;
    end
DISP1(1,t+1)= e1*SIG*(e1)';
    

end


 figure
plot(DISP1)

%normal plotter
fDISP1=zeros(200,10);
axmin=max(DISP1.^.5)*4;
for jj=1:10;
for j=1:200;
 sig=DISP1(1,jj); 
   x=-axmin+j*axmin*0.01;
 
    fy(j,1)=(1/((sig*(2*pi)^.5)))*exp(-(x^2)/(2*sig^2));  
    fDISP1(j,jj)=fy(j,1);
    
end


end
figure
surf(fDISP1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARCH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
periods=20;
ZVEC=[zeros(1, bindim-1) 1 zeros(1,periods-1) ];

for t=1:periods;
    j=binbase*ZVEC(1,t:t+bindim-1)';
    z=ZVEC(1,t:t+bindim-1);
    if z(end) ==1;
        SENS(t)=A(1,:,j+1)*N(:,1,j+1);
    else
        SENS(t)=theta(6)*A(1,:,j+1)*N(:,1,j+1);
    end
end
figure
plot(SENS)


