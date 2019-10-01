clear
clc
%% Sample the function f(r,theta)=h*g
N2=41 ; %number of sample points in angular direction
N1=383; %number of sample points in radial direction
M=(N2-1)/2; %highest order of bessel function
R=40;% space limit
Wp=15; % band limit
a=0.1;
load('zeromatrix.mat')
theta=thetamatrix_SpaceLimited(N2,N1);
r=rmatrix_SpaceLimited(N2,N1,R,zeromatrix);
psi=psimatrix_SpaceLimited(N2,N1);
rho=rhomatrix_SpaceLimited(N2,N1,R,zeromatrix);
[x,y]=pol2cart(theta,r); %sample points in Cartesian coordinates
[x1,y1]=pol2cart(psi,rho);

%creating a discrete true function
for ii=1:N2
    for jj=1:N1-1
        trueFunc(ii,jj)=(-6*pi*j*sin(psi(ii,jj))*(((sqrt(rho(ii,jj)^2+a^2)-a)^1)/(rho(ii,jj)^1*sqrt(rho(ii,jj)^2+a^2))))+(2*pi*j*sin(3*psi(ii,jj))*(((sqrt(rho(ii,jj)^2+a^2)-a)^3)/(rho(ii,jj)^3*sqrt(rho(ii,jj)^2+a^2))))+(-8*pi*cos(10*psi(ii,jj))*(((sqrt(rho(ii,jj)^2+a^2)-a)^10)/(rho(ii,jj)^10*sqrt(rho(ii,jj)^2+a^2))))+(24*pi*j*sin(15*psi(ii,jj))*(((sqrt(rho(ii,jj)^2+a^2)-a)^15)/(rho(ii,jj)^15*sqrt(rho(ii,jj)^2+a^2))));
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
        f(ii,jj)=(exp(-a*r(ii,jj))/r(ii,jj))*(3*sin(theta(ii,jj))+sin(3*theta(ii,jj))+4*cos(10*theta(ii,jj))+12*sin(15*theta(ii,jj)));
    end
end

error= 20*log10(abs(f- TwoDIFT)/max(max(abs(TwoDIFT))));

figure(1)
subplot(2,1,1)
surf(x,y,abs(f))
set(gca,'linewidth',1,'fontsize',25,'fontname','Times');
title('\fontsize{24}Continuous Inverse Transform')
subplot(2,1,2)
surf(x,y,abs(TwoDIFT))
set(gca,'linewidth',1,'fontsize',25,'fontname','Times');
title('\fontsize{24}Discrete inverse Transform')

figure(2)
surf(x,y,error)
set(gca,'linewidth',1,'fontsize',25,'fontname','Times');
xlabel('x');
ylabel('y');
zlabel('db')
str=sprintf('Error distribution with N2 = %d, N1 = %d,R= %d, a= %d ', N2,N1,R,a);
title(['\fontsize{24}Error distribution with N2=',num2str(N2),', N1=',num2str(N1),', R=',num2str(R), ', Wp=',num2str(Wp)]);

mean=mean(mean(error));
max=max(max(error));