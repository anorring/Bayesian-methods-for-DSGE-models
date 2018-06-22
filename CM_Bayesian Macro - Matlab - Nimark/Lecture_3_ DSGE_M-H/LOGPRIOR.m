function LP=LOGPRIOR(theta)


LP=0;
LP=LP+lpdfNormal(theta(2),2,0.1);

LP=LP+lpdfBeta(theta(3),0.75,0.05);

LP=LP+lpdfBeta(theta(4),0.99,0.1);