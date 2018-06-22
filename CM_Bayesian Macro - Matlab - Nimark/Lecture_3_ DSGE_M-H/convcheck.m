function x=convcheck(MCMC)
n=size(MCMC,1);
m=size(MCMC,2);
x=[];
for j=1:10000:m
    X=diag(cov(MCMC(:,100:j)'));
    
    x=[x X];
end
%
sqrn=n^.5;
figure
for j=1:n;
    subplot(ceil(sqrn),ceil(sqrn),j);
    plot(x(j,:));
end

figure
for j=1:n;
    subplot(ceil(sqrn),ceil(sqrn),j);
   plot(MCMC(j,10:end));
end