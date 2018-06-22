function x=convcheck(MCMC)
n=size(MCMC,1);
m=size(MCMC,2);
x=[];
mu=[];
for j=10:10:m
    X=diag(cov(MCMC(:,1:j)'));    
    x=[x X];
    MU=mean(MCMC(:,1:j),2);
    mu=[mu MU;];
end


sqrn=n^.5;


figure
for j=1:n;
    subplot(ceil(sqrn),ceil(sqrn),j);
    plot(MCMC(j,:),'linewidth',2);
    if j==1;
        xlabel({'c_{11} '});
    end;
    if j==2;
        xlabel('c_{21} ');
    end;
    if j==3;
        xlabel('c_{22} ');
    end;
    if j==4;
        xlabel('a_{11} ');
    end;
    if j==5;
        xlabel('a_{12} ');
    end;
    if j==6;
        xlabel('a_{21} ');
    end;
    if j==7;
        xlabel('a_{22}');
    end;
  

end
figure
for j=1:n;
    subplot(ceil(sqrn),ceil(sqrn),j);
    plot(mu(j,:),'linewidth',2);
    if j==1;
        xlabel({'c_{11} '});
    end;
    if j==2;
        xlabel('c_{21} ');
    end;
    if j==3;
        xlabel('c_{22} ');
    end;
    if j==4;
        xlabel('a_{11} ');
    end;
    if j==5;
        xlabel('a_{12} ');
    end;
    if j==6;
        xlabel('a_{21} ');
    end;
    if j==7;
        xlabel('a_{22}');
    end;
    title('Recursive mean of MCMC')
end
figure
for j=1:n;
    subplot(ceil(sqrn),ceil(sqrn),j);
    plot(x(j,:),'linewidth',2);
    if j==1;
        xlabel({'c_{11} '});
    end;
    if j==2;
        xlabel('c_{21} ');
    end;
    if j==3;
        xlabel('c_{22} ');
    end;
    if j==4;
        xlabel('a_{11} ');
    end;
    if j==5;
        xlabel('a_{12} ');
    end;
    if j==6;
        xlabel('a_{21} ');
    end;
    if j==7;
        xlabel('a_{22}');
    end;
    title('Recursive variance of MCMC')

end

