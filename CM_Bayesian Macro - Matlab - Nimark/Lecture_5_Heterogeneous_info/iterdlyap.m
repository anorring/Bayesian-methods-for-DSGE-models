function [S]=iterdlyap(A,CC);
S=zeros(size(A));
diff=1;tol=0.00000001;
while diff > tol;
    Sst=CC+A*S*A';
    diff=max(max(abs(Sst-S)));
    S=.5*Sst+.5*S;
end;
    