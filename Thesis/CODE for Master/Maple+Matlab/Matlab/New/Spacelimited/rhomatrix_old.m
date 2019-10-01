function rho=rhomatrix_old(N2,N1,R)
rho=zeros(N2,N1-1);
for ii=1:N2;
    q=ii-1;  
    for l=1:N1-1;
        zero2=besselzero(q,N1,1);   
        jql=zero2(l);
        rho(ii,l)=jql/R;
    end
end