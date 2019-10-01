% N1 sample size in radial direction
% N2 sample size in angular direction
function theta=thetamatrix_SpaceLimited(N2,N1)
theta=zeros(N2,N1-1);
M=(N2-1)/2;
for ii=1:N2;
    p=ii-1-M;  
    for k=1:N1-1;
        theta(ii,k)=(p/N2)*2*pi;
    end
end