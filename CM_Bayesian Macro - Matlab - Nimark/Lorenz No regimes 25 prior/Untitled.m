
ST=zeros(116,1);
for j=1:length(resid)
 if abs(resid(j)) >= 1.5*(cov(resid))^.5
     ST(j)=1;
 end
end
    
   subplot(2,1,1);plot(ST) 
   subplot(2,1,2);plot(resid) 