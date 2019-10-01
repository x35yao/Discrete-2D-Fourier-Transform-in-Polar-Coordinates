function r=rmatrix_old(N2,N1,R)
r=zeros(N2,N1-1);
for ii=1:N2;
    p=ii-1;  
    for k=1:N1-1;
        zero2=besselzero(p,N1,1);   
        jpk=zero2(k);
        jpN1=zero2(N1);
        r(ii,k)=(jpk/jpN1)*R;
    end
end
