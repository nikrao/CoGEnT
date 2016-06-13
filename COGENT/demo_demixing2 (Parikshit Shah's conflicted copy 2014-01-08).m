% test CoGEnT for demixing astronomical images
clear;
clc;
close all;
warning off
dct = false;
slr = true;

% SPARSE + LOW RANK
n=20;
r=2;
s=10;
A=randn(n,n);
matsize=size(A);
[U, S, V]=svd(A);
x2_true=U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
B=randn(n,n);
tmp1=sort(abs(B(:)),'descend');
tmp2=tmp1(s);
x1_true=B.*(abs(B)>=tmp2);
y=x1_true+x2_true;
Phi=eye(n^2,n^2);
tau2=trace(S(1:r,1:r));
tau1=norm(x1_true(:),1);
tmp3=randn(n,1);
tmp3=tmp3/norm(tmp3);
Ainit2=tmp3*tmp3';
Ainit1=zeros(n,n);
Ainit1(1,1)=1;

Ainit1=Ainit1(:);
Ainit2=Ainit2(:);
y=y(:);
selfun1 = @(gradf) find_next_atom_l1(gradf);
selfun2 = @(gradf) find_next_atom_nucnorm(gradf,matsize(1),matsize(2));

[x1,x2,At1, At2 iter, obj, time, back_count] = CoGEnT_Demix_lowrank2(y, Phi, tau1, tau2,...
    Ainit1, Ainit2, selfun1, [matsize],...
    'tol',1e-4,...
    'verbose',0,...
    'maxiter',1000,...
    'dropcount',100);
norm(x1-x1_true(:))
norm(x2-x2_true(:))
            