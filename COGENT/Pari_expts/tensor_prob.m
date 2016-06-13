clear
clc

%using tensor toolbox. See for notes:
%http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.103.8658&rep=rep1&type=pdf

n=10;

%Define tensor
A=randn(n,n);
[U1,S1,V1]=svd(A);
S2=diag(S1);
S=S2(1:3);
U=U1(:,1:3);

T=ktensor(S,U,U,U);

Tfull=full(T);
innerprod(Tfull,Tfull)

T2=Tfull(:);
T3=reshape(T2,10,10,10);
norm(T3-Tfull); 

% repeat to find eigenvector
[lambda_est, Uest, maxval, v]=tensor_eig(Tfull);
