function lprior=logpriorMBD(bcan)
lprior=0;
%prod persistence
pmean=0.9;
stdd=0.02;

a = (1-pmean)*pmean^2/stdd^2 - pmean;
b = a*(1/pmean - 1);

lprior= lprior + log(beta_pdf(bcan(1), a, b));

%demand persistence
pmean=0.7;
stdd=0.1;

a = (1-pmean)*pmean^2/stdd^2 - pmean;
b = a*(1/pmean - 1);

lprior= lprior + log(beta_pdf(bcan(2), a, b));



%Disutility of labor curvature
pmean=1;
stdd=0.1;
lprior= lprior + log(norm_pdf(bcan(11),pmean ,stdd^2));


%Demand elasticity
pmean=1;
stdd=0.2;
lprior= lprior + log(norm_pdf(bcan(12),pmean ,stdd^2));



%Taylor rule coeff lagged interest rate
pmean=0.5;
stdd=0.1;
lprior= lprior + log(beta_pdf(bcan(13), pmean, stdd^2));

%Taylor rule coeff on infl
pmean=1.5;
stdd=0.05;
lprior= lprior + log(norm_pdf(bcan(14), pmean, stdd^2));

%Taylor rule coeff on output
pmean=0.5;
stdd=0.1;
lprior= lprior + log(norm_pdf(bcan(15),pmean ,stdd^2));


%Calvo sticky param
pmean=0.7;
stdd=0.05;
a = (1-pmean)*pmean^2/stdd^2 - pmean;
b = a*(1/pmean - 1);
lprior= lprior + log(beta_pdf(bcan(16), a, b));

%Discount rate
pmean=0.99;
stdd=0.01;
a = (1-pmean)*pmean^2/stdd^2 - pmean;
b = a*(1/pmean - 1);
lprior= lprior + log(beta_pdf(bcan(17), a, b));
% 
% %omega (unconditional frequency)
% pmean=0.2;
% stdd=0.05;
% a = (1-pmean)*pmean^2/stdd^2 - pmean;
% b = a*(1/pmean - 1);
% lprior= lprior + log(beta_pdf(bcan(18), a, b));
