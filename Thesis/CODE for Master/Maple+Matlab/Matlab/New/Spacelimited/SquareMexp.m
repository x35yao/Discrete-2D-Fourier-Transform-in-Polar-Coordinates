clear
clc
%% Sample the function f(r,theta)=h*g
N2=61; %number of sample points in angular direction
N1=478; %number of sample points in radial direction
M=(N2-1)/2; %highest order of bessel function
R=30;% space limit
Wp=50; % band limit
a=1;
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
        f(ii,jj)=(exp(-a*r(ii,jj))/r(ii,jj))*(heaviside(theta(ii,jj)+pi/2)-heaviside(theta(ii,jj)-pi/2));
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
        trueFunc(ii,jj)=((-2*i/13)*((sqrt(rho(ii,jj)^2+1)-1)^13)*(exp(13*i*psi(ii,jj))))/((rho(ii,jj)^13)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/13)*((sqrt(rho(ii,jj)^2+1)-1)^13)*(exp(-13*i*psi(ii,jj))))/((rho(ii,jj)^13)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/11)*((sqrt(rho(ii,jj)^2+1)-1)^11)*(exp(11*i*psi(ii,jj))))/((rho(ii,jj)^11)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/11)*((sqrt(rho(ii,jj)^2+1)-1)^11)*(exp(-11*i*psi(ii,jj))))/((rho(ii,jj)^11)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/3)*((sqrt(rho(ii,jj)^2+1)-1)^3)*(exp(-3*i*psi(ii,jj))))/((rho(ii,jj)^3)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/3)*((sqrt(rho(ii,jj)^2+1)-1)^3)*(exp(3*i*psi(ii,jj))))/((rho(ii,jj)^3)*(sqrt(rho(ii,jj)^2+1)))+((-2*i)*((sqrt(rho(ii,jj)^2+1)-1))*(exp(1*i*psi(ii,jj))))/((rho(ii,jj))*(sqrt(rho(ii,jj)^2+1)))+((-2*i)*((sqrt(rho(ii,jj)^2+1)-1))*(exp(-1*i*psi(ii,jj))))/((rho(ii,jj))*(sqrt(rho(ii,jj)^2+1)))+((-2*i/5)*((sqrt(rho(ii,jj)^2+1)-1)^5)*(exp(-5*i*psi(ii,jj))))/((rho(ii,jj)^5)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/5)*((sqrt(rho(ii,jj)^2+1)-1)^5)*(exp(5*i*psi(ii,jj))))/((rho(ii,jj)^5)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/15)*((sqrt(rho(ii,jj)^2+1)-1)^15)*(exp(-15*i*psi(ii,jj))))/((rho(ii,jj)^15)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/15)*((sqrt(rho(ii,jj)^2+1)-1)^15)*(exp(15*i*psi(ii,jj))))/((rho(ii,jj)^15)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/7)*((sqrt(rho(ii,jj)^2+1)-1)^7)*(exp(7*i*psi(ii,jj))))/((rho(ii,jj)^7)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/7)*((sqrt(rho(ii,jj)^2+1)-1)^7)*(exp(-7*i*psi(ii,jj))))/((rho(ii,jj)^7)*(sqrt(rho(ii,jj)^2+1)))+pi/(sqrt(rho(ii,jj)^2+1))+((-2*i/17)*((sqrt(rho(ii,jj)^2+1)-1)^17)*(exp(-17*i*psi(ii,jj))))/((rho(ii,jj)^17)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/17)*((sqrt(rho(ii,jj)^2+1)-1)^17)*(exp(17*i*psi(ii,jj))))/((rho(ii,jj)^17)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/9)*((sqrt(rho(ii,jj)^2+1)-1)^9)*(exp(9*i*psi(ii,jj))))/((rho(ii,jj)^9)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/9)*((sqrt(rho(ii,jj)^2+1)-1)^9)*(exp(-9*i*psi(ii,jj))))/((rho(ii,jj)^9)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/19)*((sqrt(rho(ii,jj)^2+1)-1)^19)*(exp(19*i*psi(ii,jj))))/((rho(ii,jj)^19)*(sqrt(rho(ii,jj)^2+1)))+((-2*i/19)*((sqrt(rho(ii,jj)^2+1)-1)^19)*(exp(-19*i*psi(ii,jj))))/((rho(ii,jj)^19)*(sqrt(rho(ii,jj)^2+1)));

    end


end

%calculating the error from transform and true function using dynamic error
error= 20*log10(abs(trueFunc- TwoDFT)/max(max(abs(TwoDFT))));

figure(1)
subplot(2,1,1)
surf(x1,y1,abs(trueFunc))
title('\fontsize{24}Sampled Continuous Forward Transform')
subplot(2,1,2)
surf(x1,y1,abs(TwoDFT))
title('\fontsize{24}Discrete Forward Transform')

figure(2)

surf(x1,y1,error)
xlabel('x');
ylabel('y');
zlabel('db')
str=sprintf('Error distribution with N2 = %d, N1 = %d,R= %d, a= %d ', N2,N1,R,a);
title(['\fontsize{24}Error distribution with N2=',num2str(N2),', N1=',num2str(N1),', R=',num2str(R), ', Wp=',num2str(Wp)]);

mean1=mean(mean(error));
max1=max(max(error));
