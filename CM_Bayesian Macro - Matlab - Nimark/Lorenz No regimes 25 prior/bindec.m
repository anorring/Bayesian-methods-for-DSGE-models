%binary vector - decimal conversions


%decimal to binary vector
binmaxdim=8;
n=7
v=dec2bin(n)=='1'
binvec=[zeros(1,binmaxdim-length(v)) v ;]

%binary vector to decimal

for j=1:binmaxdim
    binbase(1,binmaxdim-j+1)=2^(j-1);
end
binbase
binbase*binvec'

