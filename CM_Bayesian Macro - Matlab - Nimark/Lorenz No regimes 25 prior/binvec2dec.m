function d=binvec2dec(v)

for j=1:length(v)
    binbase(1,length(v)-j+1)=2^(j-1);
end
d=binbase*v';
