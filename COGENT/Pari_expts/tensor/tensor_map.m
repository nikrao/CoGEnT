function [X, w, lambda ] = tensor_map( Tfull,v )
%Computes X=T(I,I,v) w=T(I,v,v) lambda=T(v,v,v) as defined in Anandkumar et al.
%Here T is a symmetric tensor of order 3 and dimension n.
n=length(v);

w=zeros(n,1);
I=eye(n,n);
%for i=1:n
%    for j=1:n
%        for l=1:n
%            w=w+Tfull(i,j,l)*v(j)*v(l)*I(:,i);
%        end
%    end
%end
X=ttv(Tfull,v,3);
w=ttv(X,v,2);
lambda=ttv(w,v,1);
X=double(X);
w=double(w);
lambda=double(lambda);


end

