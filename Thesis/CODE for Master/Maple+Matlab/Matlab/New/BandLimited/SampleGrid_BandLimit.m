clc;
N1=96;
N2=95;
Wp=90;
load('zeromatrix.mat')
%% sampling grids for regular space 
theta=thetamatrix_BandLimited(N2,N1);
r=rmatrix_BandLimited(N2,N1,Wp,zeromatrix);
figure(1)
a=polar(theta,r,'.r')

%% sampling grids for frequency space
psi=psimatrix_BandLimited(N2,N1);
rho=rhomatrix_BandLimited(N2,N1,Wp,zeromatrix);
figure(2)
a=polar(psi,rho,'.r')
 
