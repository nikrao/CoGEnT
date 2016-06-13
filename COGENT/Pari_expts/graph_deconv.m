%Graph deconvolution
clear
clc
n=50;

%Generate the tree
Tree=rand(n,n);
Tree=Tree+Tree';
[T,C]=UndirectedMaximumSpanningTree(Tree);
A1=T.*Tree;
A1=A1/norm(A1,'fro');
A1=2*A1;
%A1=zeros(n,n);
%view(biograph(A1))

%Generate the cycle like graph
I=eye(n,n);
P=circulant(I(2,:));
C=rand(n,n).*P;
C(floor(n/2),n)=rand(1,1);
A2=C+C';
view(biograph(A2))
%A2=zeros(n,n);

%Their sum
A=A1+A2;
%view(biograph(A))

%Deconvolution
%Atomic set 1 = Trees
%Atomis set 2 = permutations of chorded cycle graph

% COGENT Formulation
N=length(A(:));
Phi=eye(N,N);
tau1=norm(A1,'fro');
tau2=norm(A2,'fro'); 
%tau2=.0001;
Ainit1=UndirectedMaximumSpanningTree(rand(n,n));
Ainit1=Ainit1/norm(Ainit1,'fro');
Ainit2=(C>0);
Ainit2=P^2*Ainit2*(P^2)';
Ainit2=Ainit2+Ainit2';
Ainit2=Ainit2/norm(Ainit2,'fro');
y=A(:);
Ainit1=Ainit1(:);
Ainit2=Ainit2(:);

selfun1 = @(gradf) selfun_tree(gradf);
selfun2 = @(gradf) selfun_circ(gradf);
[x1,x2,At1, At2 iter, obj, time, back_count]=CoGEnT_Demix2(y, Phi, tau1, tau2, Ainit1, Ainit2, selfun1, selfun2,size(A1),'maxiter',200);
A1_hat=real(reshape(x1,n,n));
A1_hat=A1_hat.*(A1_hat>=.09);
A2_hat=real(reshape(x2,n,n));
A2_hat=A2_hat.*(A2_hat>=.09);
view(biograph(A1))
view(biograph(A2))
view(biograph(A1_hat))
view(biograph(A2_hat))


