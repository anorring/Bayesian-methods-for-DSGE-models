function x=convcheck(MCMC)
n=size(MCMC,1);
m=size(MCMC,2);
x=[];
mu=[];
for j=1:10000:m
    X=diag(cov(MCMC(:,10:j)'));    
    x=[x X];
    MU=mean(MCMC(:,1:j),2);
    mu=[mu MU];
end
%
sqrn=n^.5;
figure(1)
for j=1:n;
    subplot(ceil(sqrn),ceil(sqrn),j);
    plot(x(j,:));
end

figure(2)
for j=1:n;
    subplot(ceil(sqrn),ceil(sqrn),j);
   plot(MCMC(j,10:end));
end

figure(3)
for j=1:n;
    subplot(ceil(sqrn),ceil(sqrn),j);
   plot(mu(j,:));
end