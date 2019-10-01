clear
clc
%% Forward Tranform
N2=15; %number of sample points in angular direction
N1=383; %number of sample points in radial direction
M=(N2-1)/2; %highest order of bessel function
R=40;% space limit
Wp=30; % band limit
a=0.1;
load('zeromatrix.mat')
theta=thetamatrix_SpaceLimited(N2,N1);
r=rmatrix_SpaceLimited(N2,N1,R,zeromatrix);
psi=psimatrix_SpaceLimited(N2,N1);
rho=rhomatrix_SpaceLimited(N2,N1,R,zeromatrix);
[x,y]=pol2cart(theta,r); %sample points in Cartesian coordinates
[x1,y1]=pol2cart(psi,rho);

%dicretizing the function
for ii=1:N2
    for jj=1:N1-1
        f(ii,jj)=exp(-r(ii,jj)^2);
    end
end

fnk=circshift(fft(circshift(f,M+1,1),N2,1),-(M+1),1);

for n=-M:M
    ii=n+M+1;
    zero2=zeromatrix(5001-abs(n),:);
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
%% Inverse Transform
FNL=circshift(fft(circshift(TwoDFT,M+1,1),N2,1),-(M+1),1);

for n=-M:M
    ii=n+M+1;
    zero2=zeromatrix(5001-abs(n),:);
    jnN1=zero2(N1);
    if n<0
    Y=((-1)^abs(n))*YmatrixAssembly(abs(n),N1,zero2);
    else
    Y=YmatrixAssembly(abs(n),N1,zero2);
    end
    Y0=Y';
    Fnk(ii,:)=FNL(ii,:)*Y0;
    fnk(ii,:)=Fnk(ii,:)*((jnN1)*(j^n))/(2*pi*R^2);
    
end
TwoDIFT=circshift(ifft(circshift(fnk,M+1,1),N2,1),-(M+1),1);
%% Precision
precision=sum(sum(abs(f-TwoDIFT)))/(N1*N2)
