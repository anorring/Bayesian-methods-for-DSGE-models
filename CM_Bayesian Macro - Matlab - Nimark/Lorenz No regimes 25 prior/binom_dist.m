function pmf = binom_dist(N,p,k)
  nValues = numel(k);
  pmf = zeros(1,nValues);
  for i = 1:nValues
    pmf(i) = nchoosek(N,k(i))*p^k(i)*(1-p)^(N-k(i));
  end
end