% Three ways to solve a model under full information rational expectations (c) Kristoffer Nimark
%K Nimark, Barcelona GSE Summer School 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign values to parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=.99; %  Discount rate
f=1.5 ; % Taylor rule paramter
k=.5; %   Slope of Phillips Curve
r=.9; %   Persistence of potential output
s=2;%     Standard deviation of shocks to potential output


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method #1
% Soderlind of stable/unstable decoupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic %start timer
cutoff=.999999; %Define the cutoff for stable vs unstable eignevalues (should be just below unity).

% Define model matrices
A0=[1,0,0;0,b,0;0,s,1;];
A1=[r,0,0;k,1,-k;0,s*f,1;];
A=A0\A1;

egen = abs(eig(A)) < cutoff; %%% defining the eigenvalues: does the model have a solution?

n1=1; %Number of predetermined variables
n2=2; % Number of jump variables
n = n1 + n2; %Total number of variables

%MatLab, complex generalized Schur decomposition
%%% matlab cannot do the analytical Schur decomposition, but there is a
%%% work-a-around using qz:
[S,T,Qa,Z] = qz(eye(size(A)),A); %MatLab: I=Q'SZ' and A=Q'TZ'; Paul S:  I=QSZ' and A=QTZ',%but Q isn't used
[S,T,Qa,Z] = reorder(S,T,Qa,Z);   % reordering of generalized eigenvalues, T(i,i)/S(i,i), in ascending order
logcon = abs(diag(T)) <= (abs(diag(S))*cutoff);  %1 for stable eigenvalue

if sum(logcon) < n1;
    warning('Too few stable roots: no stable solution');
    M = NaN; C = NaN;J0 = NaN;
    return;
elseif sum(logcon) > n1;
    warning('Too many stable roots: inifite number of stable solutions');
    M = NaN; C = NaN;J0 = NaN;
    return;
end;


Stt = S(1:n1,1:n1);
Zkt = Z(1:n1,1:n1);
Zlt = Z(n1+1:n,1:n1);
Ttt = T(1:n1,1:n1);

if cond(Zkt) > 1e+14; 
    warning('Zkt is singular: rank condition for solution not satisfied');
    M = NaN; C = NaN;J0 = NaN;
    return;
end;

Zkt_1 = inv(Zkt);         %inverting
Stt_1 = inv(Stt);


M = real(Zkt*Stt_1*Ttt*Zkt_1);       %x1(t+1) = M*x1(t) + e(t+1), %%% should be equal to rho defined above
G = real(Zlt*Zkt_1) ;              %x2(t) = C*x1(t), %%% a and b in the equation for inflation and output gap


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display solution/output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc %display time passed since "tic"

display('Stable/unstable eigenvalue decoupling');

G


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method #2
% Undetermined coefficents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
deno=r+b*r-b*r*r-k*s*f+k*s*r-1; %%% deno = denominator
a=-k*(r-1)/deno;
bb=-k*s*(f-r)/deno;
toc
display('Undetermined coefficent solution [a bb;]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display solution/output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a, bb;]'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method #3
% Project expected inflation and output on current inflation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

%starting values
c0=1;
d0=-.3;
%keep iterations (for plotting) in these vectors 
resc=[];
resd=[];

DIFF=1;
%start the  main loop
while DIFF > 0.00000001;
    A0=[1,0,0;k,1-b*c0,-k;0,-d0+s*f-s*c0,1];
    A1=[r,0,0;0,0,0;0,0,0;];
    A=A0\A1;
    C=A0\[1,0,0;]';


    S=dlyap(A,C*C'); %S is unconditional covariance of state 
    AS=A*S;

    c1=AS(2,2)/S(2,2);
    d1=AS(2,3)/S(2,2);

    DIFF=max(abs([c0 d0]-[c1 d1])); %%% c0 is a initial guess
    resc=[resc c0];
    resd=[resd d0];

    c0=c1;
    d0=d1;
end;

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display solution/output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Iterated linear projections solution [c d]');

C(2:3)



% Uncomment these lines to plot iterations of c and d. 
% Change paramter r (exogenous persistence) and see how it affects
% convergence rate
figure
subplot(2,1,1);plot(resc);
%legend('c(s)','Fontsize',16);
subplot(2,1,2);plot(resd);
%legend('d(s)','Fontsize',16);
title('Speed of convergence', 'Fontsize',24)