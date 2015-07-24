function atom = find_next_atom_nucnorm(gradf,n1,n2)

%find the next atom for nuclear norm constrained problem

% gradf is the gradient "vector", corresponding to the derivative of the 
% frobenius norm wrt X

%convert gradient to a matrix
gradf = reshape(gradf,n1,n2);


[U,~,V] = lansvd(-gradf,1,'L'); % compute largest singular vectors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [U,~,V] = svds(-gradf,1);
% % for smaller/non low rank matrices, use the next line instead of the above
% [U,~,V] = svd(-gradf,'econ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = V';
atom = U(:,1)*V(1,:);
atom = atom(:)/norm(atom(:));
clear U V;

end