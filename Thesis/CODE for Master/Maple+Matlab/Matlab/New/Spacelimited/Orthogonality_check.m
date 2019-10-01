%% Othogonality Check 1
clc;
N1=16; % sampling points in radial direction
N2=15; % sampling points in angular direction
M=(N2-1)/2;
A=zeros(N2,N1-1); % define a N2*(N1-1)zero matrix
p_0=5; % p'
k_0=11; % k'
q_0=5; % q'
l_0=11; % l'

% Calculate the Bessel zeros
for n=-M:M;
    a=n+M+1;
        C=besselzero(n,N1,1); 
        zero(a,:)=C;    
end
% 
for l=1:N1-1;
    for ii=1:N2;
        q=ii-1-M;
        A1=kernelminus(N2,N1,q,l,zero); % kernelmius with entry p,k
        A2=kernelplus(N2,N1,q,l,zero); % kernelplus with entry p,k
        A=A+A1*A2(p_0,k_0);
       
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %Othogonality Check 2
B=zeros(N2,N1-1);% Define a N2*(N1-1)zero matrix

for k=1:N1-1;
    for ii=1:N2;
        p=ii-1-M;
        B1=kernelminus2(N2,N1,p,k,zero);% kernelmius with entry q,l
        B2=kernelplus2(N2,N1,p,k,zero); %kernelplus with entry q,l
        B=B+B1*B2(q_0,l_0);
       
    end
end

      
             
        
        