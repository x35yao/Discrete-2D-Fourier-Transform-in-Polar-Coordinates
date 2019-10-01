clear
clc
%% Sample the function f(r,theta)=h*g
N1_vector=[190,240,290,340,390];
N2_vector=[7,15,31,61,121];
E_max=zeros(5,5);
E_avg=zeros(5,5);
for mm=1:5
    for nn=1:5   
N2=N2_vector(mm); %number of sample points in angular direction
N1=N1_vector(nn); %number of sample points in radial direction
M=(N2-1)/2; %highest order of bessel function
R=150;% space limit
Wp=6; % band limit
a=0.1;
load('zeromatrix.mat')
theta=thetamatrix_SpaceLimited(N2,N1);
r=rmatrix_SpaceLimited(N2,N1,R,zeromatrix);
psi=psimatrix_SpaceLimited(N2,N1);
rho=rhomatrix_SpaceLimited(N2,N1,R,zeromatrix);
[x,y]=pol2cart(theta,r); %sample points in Cartesian coordinates
[x1,y1]=pol2cart(psi,rho);

clear f
clear trueFunc
fnl=zeros(N2,N1-1);
Fnl=zeros(N2,N1-1);
error=zeros(N2,N1-1);

%dicretizing the function
for ii=1:N2
    for jj=1:N1-1
        f(ii,jj)=heaviside(10-r(ii,jj))-heaviside(5-r(ii,jj));
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
    Y0=Y';
    fnl(ii,:)=fnk(ii,:)*Y0;
    Fnl(ii,:)=fnl(ii,:)*(2*pi*(i^(-n)))*(R^2/jnN1);
end
TwoDFT=circshift(ifft(circshift(Fnl,M+1,1),N2,1),-(M+1),1);
%creating a discrete true function
for ii=1:N2
    for jj=1:N1-1
        trueFunc(ii,jj)=2*pi*((10/rho(ii,jj))*besselj(1,10*rho(ii,jj))-(5/rho(ii,jj))*besselj(1,5*rho(ii,jj)));
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
