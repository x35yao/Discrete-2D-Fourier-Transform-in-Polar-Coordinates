%Build up a matrix from 1th to 10000th zero of 
%Bessel funtions from order 0 to order 5000
N1=10000;
N2=10001;
M=(N2-1)/2;
zeromatrix=zeros(M+1,N1);
for n=-M:0;
    ii=n+M+1;
    zero1=besselzero(n,N1,1);
    zeromatrix(ii,:)=zero1;
end;
