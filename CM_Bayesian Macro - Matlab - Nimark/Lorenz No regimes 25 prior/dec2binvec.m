function v=dec2binvec(d,dim)

vv=dec2bin(d)=='1';

if d <= 2^(dim-1)
    v=[zeros(1,dim-length(vv)) vv ;];
else
    v= vv ;
end