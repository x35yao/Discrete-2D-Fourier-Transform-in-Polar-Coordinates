clear
clc
%% Sample the function f(r,theta)=h*g
N1_vector=[330,380,430,480,530];
N2_vector=[11,21,41,81,161];
E_max=zeros(5,5);
E_avg=zeros(5,5);
for mm=1:5
    for nn=1:5
N2=N2_vector(mm); %number of sample points in angular direction
N1=N1_vector(nn); %number of sample points in radial direction
M=(N2-1)/2; %highest order of bessel function
R=15;% space limit
Wp=90; % band limit
a=5;
load('zeromatrix.mat')
theta=thetamatrix_BandLimited(N2,N1);
r=rmatrix_BandLimited(N2,N1,Wp,zeromatrix);
psi=psimatrix_BandLimited(N2,N1);
rho=rhomatrix_BandLimited(N2,N1,Wp,zeromatrix);
[x,y]=pol2cart(theta,r); %sample points in Cartesian coordinates
[x1,y1]=pol2cart(psi,rho);
clear f
clear trueFunc
Fnk=zeros(N2,N1-1);
fnk=zeros(N2,N1-1);
error=zeros(N2,N1-1);

%creating a discrete true function
for ii=1:N2
    for jj=1:N1-1
        
        if rho(ii,jj)<a
        trueFunc(ii,jj)=(8*pi*cos(10*psi(ii,jj))*(rho(ii,jj)^10))/(a*sqrt(a^2-rho(ii,jj)^2)*(a+sqrt(a^2-rho(ii,jj)^2))^10);
        elseif rho(ii,jj)>a
        trueFunc(ii,jj)=(-6*a*pi*j*sin(psi(ii,jj))+2*pi*j*sin(3*asin(a/rho(ii,jj)))*rho(ii,jj)*sin(3*psi(ii,jj))-8*pi*sin(10*asin(a/rho(ii,jj)))*rho(ii,jj)*cos(10*psi(ii,jj))+24*pi*j*sin(15*asin(a/rho(ii,jj)))*rho(ii,jj)*sin(15*psi(ii,jj)))/(a*rho(ii,jj)*sqrt(rho(ii,jj)^2-a^2));
        end
    end
end
FNL=circshift(fft(circshift(trueFunc,M+1,1),N2,1),-(M+1),1);

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
    fnk(ii,:)=Fnk(ii,:)*((Wp^2)*(j^n))/(2*pi*jnN1);
    
end
TwoDIFT=circshift(ifft(circshift(fnk,M+1,1),N2,1),-(M+1),1);
%%dicretizing the function in space domain
for ii=1:N2
    for jj=1:N1-1
        f(ii,jj)=(sin(a*r(ii,jj))/(a*r(ii,jj)))*(3*sin(theta(ii,jj))+sin(3*theta(ii,jj))+4*cos(10*theta(ii,jj))+12*sin(15*theta(ii,jj)));;
    end
end
error= 20*log10(abs(f- TwoDIFT)/max(max(abs(TwoDIFT))));
max1=max(max(error));
mean1=mean(mean(error));
E_max(mm,nn)=max1;
E_avg(mm,nn)=mean1;
    end
end

