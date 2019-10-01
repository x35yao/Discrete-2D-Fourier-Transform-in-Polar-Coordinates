%% Othogonality Check 1
clc;
N1=16; % sampling points in radial direction
N2=15; % sampling points in angular direction
M=(N2-1)/2;
A=zeros(N2,N1-1); % define a N2*(N1-1)zero matrix
p_0=10;
k_0=11;
q_0=4;
l_0=11;

% Calculate the Bessel zeros
for n=-M:M;
    a=n+M+1;
        C=besselzero(n,N1,1); 
        zero(a,:)=C;    
end
% Process a sum of N2*(N1-1) matrices
for l=1:N1-1;
    for ii=1:N2;
        q=ii-1-M;
        A1=kernelminus(N2,N1,q,l,zero);
        A2=kernelplus(N2,N1,q,l,zero);
        A=A+A1*A2(p_0,k_0);
       
    end
end
%% Othogonality Check 2
B=zeros(N2,N1-1);% Define a N2*(N1-1)zero matrix
B2=kernelplus(N2,N1,q_0,l_0,zero);
for ii=1:N2;
    q=ii-1-M;
    for l=1:N1-1;
        B1=kernelminus(N2,N1,q,l,zero);
        B(ii,l)=sum(sum(B1.*B2));
    end
end

      
             
        
        