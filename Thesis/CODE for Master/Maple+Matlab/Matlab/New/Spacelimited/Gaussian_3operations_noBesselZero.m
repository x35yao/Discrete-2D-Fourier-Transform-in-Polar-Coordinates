clc;
N1=383;
N2=15;
M=(N2-1)/2;
R=40;
Wp=30;

theta=thetamatrix_SpaceLimited(N2,N1);
r=rmatrix_old(N2,N1,R);
psi=psimatrix_SpaceLimited(N2,N1);
rho=rhomatrix_old(N2,N1,R);


for ii=1:N2
    for jj=1:N1-1
        f(ii,jj)=exp(-r(ii,jj)^2);
    end
end
fnk=circshift(fft(circshift(f,M+1,1),N2,1),-(M+1),1);

for n=-M:M
    ii=n+M+1;
    zero2=besselzero(n,N1,1);
    jnN1=zero2(N1);
    if n<0
    Y=((-1)^abs(n))*YmatrixAssembly(abs(n),N1,zero2);
    else
    Y=YmatrixAssembly(abs(n),N1,zero2);
    end
    fnl(ii,:)=(Y*fnk(ii,:)')';
    Fnl(ii,:)=fnl(ii,:)*(2*pi*(i^(-n)))*(R^2/jnN1);
end

TwoDFT=circshift(ifft(circshift(Fnl,M+1,1),N2,1),-(M+1),1);
