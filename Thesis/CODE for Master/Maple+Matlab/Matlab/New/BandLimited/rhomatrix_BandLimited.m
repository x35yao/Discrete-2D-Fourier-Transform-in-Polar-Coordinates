function rho=rhomatrix_BandLimited(N2,N1,Wp,zeromatrix)
M=(N2-1)/2;
for ii=1:N2;
    q=ii-1-M;  
    for l=1:N1-1;
        zero2=zeromatrix(5001-abs(q),:); 
        jql=zero2(l);
        jqN1=zero2(N1);
        rho(ii,l)=(jql/jqN1)*Wp;
    end
end