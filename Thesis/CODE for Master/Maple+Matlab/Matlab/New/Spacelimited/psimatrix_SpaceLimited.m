% N1 sample size in radial direction
% N2 sample size in angular direction
function psi=psimatrix_SpaceLimited(N2,N1)
psi=zeros(N2,N1-1);
M=(N2-1)/2;
for ii=1:N2;
    q=ii-1-M;  
    for l=1:N1-1;
        psi(ii,l)=(q/N2)*2*pi;
    end
end