function LP=log_prior_DSGE(theta)


LP=0;

% %(1) r - productivity persistence
% pmean=0.9; pstdd=0.05;
% a = (1-pmean)*pmean^2/pstdd^2 - pmean;
% b = a*(1/pmean - 1);
% LP = LP + lpdfBeta(theta(1),a,b); 

%(2) g - Consumption utility curvature/coeff. of relative risk aversion
pmean=3; pstdd=0.05;
LP = LP + lpdfNormal(theta(2),pmean,pstdd) ;

%(3) d - Calvo price setting parameter
pmean=0.75; pstdd=0.05;
a = (1-pmean)*pmean^2/pstdd^2 - pmean;
b = a*(1/pmean - 1);

LP = LP + lpdfBeta(theta(3),a,b); 

%(4) b - discount factor
pmean=0.99; pstdd=0.01;
a = (1-pmean)*pmean^2/pstdd^2 - pmean;
b = a*(1/pmean - 1);

LP = LP + lpdfBeta(theta(4),a,b);

% 
%(5) f - coefficient on inflation in Taylor rule
pmean=1.5; pstdd=0.1;
LP = LP + lpdfNormal(theta(5),pmean,pstdd) ;