clc;
N1=96;
N2=95;
R=1;
load('zeromatrix.mat')
%% sampling grids for regular space 
theta=thetamatrix_SpaceLimited(N2,N1);
r=rmatrix_SpaceLimited(N2,N1,R,zeromatrix);
figure(1)
a=polar(theta,r,'.r')

%% sampling grids for frequency space
psi=psimatrix_SpaceLimited(N2,N1);
rho=rhomatrix_SpaceLimited(N2,N1,R,zeromatrix);
figure(2)
a=polar(psi,rho,'.r')

