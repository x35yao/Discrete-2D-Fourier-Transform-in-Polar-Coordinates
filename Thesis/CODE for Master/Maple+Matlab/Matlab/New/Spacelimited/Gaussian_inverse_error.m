clear
clc
%% Sample the function f(r,theta)=h*g
N1_vector=[283,333,383,433,483];
N2_vector=[3,7,15,31,61];
E_max=zeros(5,5);
E_avg=zeros(5,5);
for mm=1:5
    for nn=1:5   
N2=N2_vector(mm); %number of sample points in angular direction
N1=N1_vector(nn); %number of sample points in radial direction
M=(N2-1)/2; %highest order of bessel functionR=40;% space limit
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

clear f
clear trueFunc
Fnk=zeros(N2,N1-1);
fnk=zeros(N2,N1-1);
error=zeros(N2,N1-1);

%creating a discrete true function
for ii=1:N2
    for jj=1:N1-1
       trueFunc(ii,jj)=pi*exp((-rho(ii,jj)^2)/4);
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
    fnk(ii,:)=Fnk(ii,:)*((jnN1)*(j^n))/(2*pi*R^2);
    
end
TwoDIFT=circshift(ifft(circshift(fnk,M+1,1),N2,1),-(M+1),1);

%%dicretizing the function in space domain
for ii=1:N2
    for jj=1:N1-1
        f(ii,jj)=exp(-r(ii,jj)^2);
    end
end

error= 20*log10(abs(f- TwoDIFT)/max(max(abs(TwoDIFT))));
max1=max(max(error));
mean1=mean(mean(error));
E_max(mm,nn)=max1;
E_avg(mm,nn)=mean1;
    end
end
E_max_round=round(E_max,1);
E_avg_round=round(E_avg,1);
