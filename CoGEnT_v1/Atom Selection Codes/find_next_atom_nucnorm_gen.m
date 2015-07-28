function atom = find_next_atom_nucnorm_gen(gradf,n1,n2,A,B)

%find the next atom for nuclear norm constrained problem

% gradf is the gradient "vector", corresponding to the derivative of the 
% frobenius norm wrt X. A and B are matrices such that we compute the generalized SVD of AgradB
	
	

%convert gradient to a matrix
gradf = reshape(gradf,n1,n2);

G = A'*gradf*B;

[U,~,V] = lansvd(-G,1,'L'); % compute largest singular vectors

V = V';
atom = U(:,1)*V(1,:);
atom = atom(:)/norm(atom(:));
clear U V;

end