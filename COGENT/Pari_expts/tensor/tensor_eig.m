function [ lambda_est, Uest, maxeig, max_eigvec ] = tensor_eig( Tfull )
%Gives eigendecomposition of an orthogonal 3-tensor
%   Tfull is the tensor, lambda_est is vector of eigenvalues and Uest is the
%   vector of eigenvectors
% maxeig is the maximum eigenvalue and max_eigvec is the corresponding
% eigenvector
n=size(Tfull,1);
Uest=[];
lambda_est=[];
Nsample=20;
for i=1:Nsample
    v=randn(n,1);
    for j=1:100
        [X,w,lambda]=tensor_map(Tfull,v);
        w=w/norm(w);
        v=w;
    end
    if sum((svd([Uest v])>=1e-3))>=sum((svd(Uest)>=1e-3))+1
        Uest=[Uest v];
        lambda_est=[lambda_est lambda];
        
    end
end
[maxeig,maxind]=max(lambda_est);
max_eigvec=Uest(:,maxind);


end

