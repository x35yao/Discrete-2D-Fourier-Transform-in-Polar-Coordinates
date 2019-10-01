function r=rmatrix_BandLimited(N2,N1,Wp,zeromatrix)
M=(N2-1)/2;
for ii=1:N2;
    p=ii-1-M;  
    for k=1:N1-1;
        zero2=zeromatrix(5001-abs(p),:);   
        jpk=zero2(k);
        r(ii,k)=jpk/Wp;
    end
end