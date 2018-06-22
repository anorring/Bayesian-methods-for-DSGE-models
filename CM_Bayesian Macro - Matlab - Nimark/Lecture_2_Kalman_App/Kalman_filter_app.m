% Kalman filter applications
% Code to play around with Kalman filter for state space system defined as
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	X[t] = AX[t-1] + Cu[t]
%
%   Z[t] = DX[t] + Ru[t]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all

periods=200;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % %Scalar case
% A=0.95;
% C=1;
% D=1;
% R=.5;


%Bivariate state, scalar signal
% A=[0.99, 0.1; -.3,0.5;];
% C=eye(2);
% D=[1 1];
% R=1;

% %Bivariate state, two dimensional signal
% A=[0.99, 0.1; -.3,0.5;];
% C=eye(2);
% D=eye(2);
% R=eye(2);

% %Bivariate state, two dimensional signal
A=eye(9);
C=eye(9);
D=eye(9);
R=10*eye(9);

[K,P,Kss,Pss,Z]=kalman_filter_sim(A,C,D,R,periods);

[LL]=LL_state_space(A,C,D,R,Z)