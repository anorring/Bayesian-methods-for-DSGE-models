function plotpost(bb_,addplot)
%%
ysca=max(size(bb_));
Q=50; %HP smoothing parameter
% figure('yscaame','Calvo and Indexing','yscaumberTitle','off');
figure(5)
for j=1:length(bb_(:,1));
    hold on;
    [n,xout] = hist(bb_(j,:),100);
    n(1,1)=0;n(end,1)=0;
    xsca=(max(xout)-min(xout))/100;
    nn=hpfilter(n,Q);
    subplot(3,3,j);
    if addplot==1;
        AXX=axis;
        AX=[ min([min(xout),AXX(1,1)]) max([max(xout),AXX(1,2)]) 0 max([(max((nn*1.1))/(ysca*xsca)),AXX(1,4)])];
        plot(xout,nn/(ysca*xsca),'g','linewidth',2);axis(AX);
    else;
        plot(xout,nn/(ysca*xsca),'linewidth',2);axis([min(xout) max(xout) 0 max((nn*1.1))/(ysca*xsca)]);
    end
    
     if j==1;
        xlabel({'\rho '});
    end;
    if j==2;
        xlabel('\gamma ');
    end;
    if j==3;
        xlabel('Calvo ');
    end;
    if j==4;
        xlabel('\beta ');
    end;
    if j==5;
        xlabel('\phi');
    end;
    if j==6;
        xlabel('\sigma_x ');
    end;
    if j==7;
        xlabel('\sigma_y ');
    end;
    if j==8;
        xlabel('\sigma_p ');
    end;
    if j==9;
        xlabel('\sigma_r ');
    end;
%     if j==10;
%         xlabel('kappa_mi');
%     end;
%     if j==11;
%         xlabel('kappa_x');
%     end;
end
