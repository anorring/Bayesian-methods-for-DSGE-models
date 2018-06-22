%produce  graphs for distributions in introduction
clear all
close all

a=0.5; %prob of z=1
sig=.15; %s.d. innovation when z=0;
gamma=2; %s.d. multiple of standard deviation when z=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mixture distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%normal plotter
fy=[];
fx=[];
fr=[];

for j=1:200;
    
    x=-4*(gamma*sig)+j*0.04*(gamma*sig);
    
    fy(j)=a*(1/((sig*(2*pi)^.5)))*exp(-(x^2)/(2*sig^2))+...
        (1-a)*(1/(((gamma*sig)*(2*pi)^.5)))*exp(-(x^2)/(2*(gamma*sig)^2));
    fx(j)=x;
    
    fya(j)=(1/((sig*(2*pi)^.5)))*exp(-(x^2)/(2*sig^2));
    fyb(j)=(1/(((gamma*sig)*(2*pi)^.5)))*exp(-(x^2)/(2*(gamma*sig)^2));
    
    fr(j)=((1-a)*(1/(((gamma*sig)*(2*pi)^.5)))*exp(-(x^2)/(2*(gamma*sig)^2)))/(a*(1/((sig*(2*pi)^.5)))*exp(-(x^2)/(2*sig^2))+(1-a)*(1/(((gamma*sig)*(2*pi)^.5)))*exp(-(x^2)/(2*(gamma*sig)^2)));
%      fr(j)=((1/(((gamma*sig)*(2*pi)^.5)))*exp(-(x^2)/(2*(gamma*sig)^2)))/((1/((sig*(2*pi)^.5)))*exp(-(x^2)/(2*sig^2))+(1/(((gamma*sig)*(2*pi)^.5)))*exp(-(x^2)/(2*(gamma*sig)^2)));
end

figure
plot(fx,fy,'linewidth',3,'linestyle','-');
legend('p(x)');

figure
plot(fx,fy,'linewidth',3,'linestyle','-');
hold on;
plot(fx,fr,'linewidth',3,'linestyle','--');
hold on;
plot(fx,fyb,'linewidth',3,'linestyle',':');
legend('p(x)', 'p(S=1|x)','p(x|S=1)' );

figure
plot(fx,fy,'linewidth',2,'linestyle','-');
hold on;
plot(fx,fya,'linewidth',2,'linestyle','--');
hold on;
plot(fx,fyb,'linewidth',2,'linestyle',':');
legend('p(x)', 'p(x|S=1)', 'p(x|S=0)');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Standard Bayesian updating illustrated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%normal plotter
fy=[];
fx=[];
fr=[];
sigeps=.1;
mux=-0.2*sig/(sigeps+sig);

postsig=inv(1/sig + 1/sigeps);

for j=1:200;
    
    x=-4*(sig)+j*0.04*(sig);
    fy(j)=(1/((sig*(2*pi)^.5)))*exp(-(x^2)/(2*sig^2));
    fb(j)=(1/((postsig*(2*pi)^.5)))*exp(-((x-mux)^2)/(2*postsig^2));
        
    fx(j)=x;
    
   
end

figure
plot(fx,fy,'linewidth',3,'linestyle','-');
hold on;
plot(fx,fb,'linewidth',3,'linestyle','--');

legend('p(x)','p(x|x_{j})');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M-B-D Bayesian updating illustrated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%normal plotter
fy=[];
fx=[];
fr=[];
sigeps=1.5;
g=((gamma*sig)^2)/((gamma*sig)^2+sigeps^2);

mux=-10*(gamma*sig)*g;

postsig=inv(1/(gamma*sig)+1/(sigeps^2));

for j=1:200;
    
    x=-4*(gamma*sig)+j*0.04*(gamma*sig);
    fy(j)=a*(1/((sig*(2*pi)^.5)))*exp(-(x^2)/(2*sig^2))+...
        (1-a)*(1/(((gamma*sig)*(2*pi)^.5)))*exp(-(x^2)/(2*(gamma*sig)^2));
    fx(j)=x;
    
%     fya(j)=(1/((sig*(2*pi)^.5)))*exp(-(x^2)/(2*sig^2));
    fyb(j)=(1/(((gamma*sig)*(2*pi)^.5)))*exp(-(x^2)/(2*(gamma*sig)^2));
    fb(j)=(1/((postsig*(2*pi)^.5)))*exp(-((x-mux)^2)/(2*postsig^2));
    
   
end

figure
plot(fx,fy,'linewidth',3);
legend('p(x)');
% hold on;
% plot(fx,fyb,'linewidth',3,'linestyle','--');
% legend('p(x)','p(x|S)');
hold on
plot(fx,fb,'linewidth',3,'linestyle',':');
% legend('p(x)','p(x|S)','p(x|y,S)');
legend('p(x)','p(x|y,S)');
