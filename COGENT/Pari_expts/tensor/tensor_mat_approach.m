clear
clc
n=10;
% define random tansor
A=randn(n,n);
[U1,S1,V1]=svd(A);
S2=diag(S1);
S=S2(1:3);
U=U1(:,1:3);
V=V1(:,1:3);
W=U1(:,4:6);

T=ktensor(S,U,V,W);

a=randn(n,1);
a=a/norm(a);

b=randn(n,1);
b=b/norm(b);

T2=full(T);

% define Ta and Tb
Ta=zeros(n,n);
Tb=zeros(n,n);
for i=1:n
    Ta=Ta+a(i)*(T2(:,:,i));
    Tb=Tb+b(i)*(T2(:,:,i));
end

T1=double(Ta)*pinv(double(Tb));
T2=double(Tb)*pinv(double(Ta));

eig(T1)
eig(T2)
