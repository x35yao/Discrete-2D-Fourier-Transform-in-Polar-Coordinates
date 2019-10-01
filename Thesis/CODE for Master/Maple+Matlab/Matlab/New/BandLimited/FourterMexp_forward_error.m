clear
clc
%% Sample the function f(r,theta)=h*g
N1_vector=[283,333,383,433,483];
N2_vector=[11,21,41,81,161];
E_max=zeros(5,5);
E_avg=zeros(5,5);
for mm=1:5
    for nn=1:5
N2=N2_vector(mm); %number of sample points in angular direction
N1=N1_vector(nn); %number of sample points in radial direction
M=(N2-1)/2; %highest order of bessel function
R=40;% space limit
Wp=30; % band limit
a=0.1;
load('zeromatrix.mat')
theta=thetamatrix_BandLimited(N2,N1);
r=rmatrix_BandLimited(N2,N1,Wp,zeromatrix);
psi=psimatrix_BandLimited(N2,N1);
rho=rhomatrix_BandLimited(N2,N1,Wp,zeromatrix);
[x,y]=pol2cart(theta,r); %sample points in Cartesian coordinates
[x1,y1]=pol2cart(psi,rho);

clear f
clear trueFunc
Fnl=zeros(N2,N1-1);
fnl=zeros(N2,N1-1);
error=zeros(N2,N1-1);
%dicretizing the function
for ii=1:N2
    for jj=1:N1-1
        f(ii,jj)=(exp(-a*r(ii,jj))/r(ii,jj))*(3*sin(theta(ii,jj))+sin(3*theta(ii,jj))+4*cos(10*theta(ii,jj))+12*sin(15*theta(ii,jj)));
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
%creating a discrete true function
for ii=1:N2
    for jj=1:N1-1
        trueFunc(ii,jj)=(-6*pi*j*sin(psi(ii,jj))*(((sqrt(rho(ii,jj)^2+a^2)-a)^1)/(rho(ii,jj)^1*sqrt(rho(ii,jj)^2+a^2))))+(2*pi*j*sin(3*psi(ii,jj))*(((sqrt(rho(ii,jj)^2+a^2)-a)^3)/(rho(ii,jj)^3*sqrt(rho(ii,jj)^2+a^2))))+(-8*pi*cos(10*psi(ii,jj))*(((sqrt(rho(ii,jj)^2+a^2)-a)^10)/(rho(ii,jj)^10*sqrt(rho(ii,jj)^2+a^2))))+(24*pi*j*sin(15*psi(ii,jj))*(((sqrt(rho(ii,jj)^2+a^2)-a)^15)/(rho(ii,jj)^15*sqrt(rho(ii,jj)^2+a^2))));
    end
end

%calculating the error from transform and true function using dynamic error
error= 20*log10(abs(trueFunc- TwoDFT)/max(max(abs(TwoDFT))));
max1=max(max(error));
mean1=mean(mean(error));
E_max(mm,nn)=max1;
E_avg(mm,nn)=mean1;
    end
end