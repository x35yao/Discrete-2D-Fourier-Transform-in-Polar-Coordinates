function A=kernelminus2(N2,N1,p,k,zero)
% Create a N2*(N1-1) zero matrix.
A=zeros(N2,N1-1);
% Calculate the order of the bessel fuction.
M=(N2-1)/2;
% Calculate the value of each entry of the matrix.
for ii=1:N2
    q=ii-1-M;
    for jj=1:N1-1
        l=jj;
        for n=-M:M;
        a=n+M+1; 
        jnk=zero(a,k);
        jnl=zero(a,l);
        jnN1=zero(a,N1);
        A(ii,jj)=A(ii,jj)+((2/N2)*(besselj(n,(jnk*jnl/jnN1))/((jnN1)*(besselj(n+1,jnk)^2)))*(j^(-n))*(exp(-j*n*2*pi*p/N2))*(exp(j*n*2*pi*q/N2)));
        end
          
    end
    
     
end   