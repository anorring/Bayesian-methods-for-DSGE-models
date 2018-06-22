function [M,L,U]=dist_plot(A,C,D,R,Z,upper,lower,ndraws,plotplease,legendplease)

%INPUT
%     Xt = A*Xt-1 +C*ut
%
%     Zt=  D1*Xt + R*ut

% 'upper' and 'lower' are percentiles (e.g. 0.975 and 0.025)
%  to be plotted along with median
%ndraws are number of draws used to construct distribution
% set plotplease = 1 if you want plots
% set legendplease = 1 if you want legends

%OUTPUT
% M,L,U are matrices with rows containing the median,
% lower and upper percentiles of X

low=ceil(ndraws*lower);
median=ceil(ndraws*.5);
upp=ceil(ndraws*upper);
dimX=length(A);
T=length(Z);
Xdist=zeros(dimX,T+1,ndraws);

for s=1:ndraws
    Xdist(:,:,s)=sim_smooth(A,C,D,R,Z);
end
Xsort=sort(Xdist,3);

if plotplease==1;
    
    q=ceil(dimX^.5);
    M=reshape(Xsort(:,:,median),dimX,T+1);
    L=reshape(Xsort(:,:,low),dimX,T+1);
    U=reshape(Xsort(:,:,upp),dimX,T+1);
    figure
    for jj=1:dimX
        subplot(q,q,jj);
        plot(L(jj,:),'linewidth',2,'linestyle','--','color','k');
        hold on
        plot(M(jj,:),'linewidth',2,'linestyle','-','color','b');
        hold on
        plot(U(jj,:),'linewidth',2,'linestyle','--','color','k');
        if legendplease==1
            legend('lower','median','upper');
        end
    end
end

