clc
clear
N1=43;  %number of sample points in radial directio
N2=15; %number of sample points in angular direction
M=(N2-1)/2;
R=40;% space limit
Wp=30;% band limit
load('zeromatrix.mat')
theta=thetamatrix_SpaceLimited(N2,N1);  %Sample point in angular direction in space domain.
r=rmatrix_SpaceLimited(N2,N1,R,zeromatrix); %Sample point in radial direction in space domain.
psi=psimatrix_SpaceLimited(N2,N1); %Sample point in angular direction in frequency domain.
rho=rhomatrix_SpaceLimited(N2,N1,R,zeromatrix); %Sample point in radial direction in frequency domain.
f=zeros(N2,N1-1); %sample points in Cartesian coordinates in space domian
F=zeros(N2,N1-1);%sample points in Cartesian coordinates in frequency domian

% Calculating Bessel zeros
for n=-M:M;
    a=n+M+1;
        C=besselzero(n,N1,1); 
        zero(a,:)=C;
end

%Dicretizing the function
for ii=1:N2
    for jj=1:N1-1
        f(ii,jj)=exp(-r(ii,jj)^2);
    end
end

for ii=1:N2;    
    q=ii-M-1;
    for jj=1:N1-1;
        l=jj;
       A=kernelminus(N2,N1,q,l,zero);
       F(ii,jj)=sum(sum(A.*f));


    end
end


