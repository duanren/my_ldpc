function Signal=mybec(modSignal,epsilon)
%本函数是对BEC信道的实现
%即在BPSK调制下，被擦除的符号的幅值视为0
if epsilon<0
    error('wrong p');
end
[m,n]=size(modSignal);
Signal=modSignal;
xx=rand(m,n);
for ii=1:m
    for jj=1:n
        if xx(ii,jj)<epsilon
            Signal(ii,jj)=0;
        end
    end
end