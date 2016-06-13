function atom = find_next_atom_tensor(gradf,n)

%find the next atom for nuclear norm constrained problem

% gradf is the gradient "vector", corresponding to the derivative of the 
% frobenius norm wrt X

%convert gradient to a tensor
Tfull = reshape(gradf,n,n,n);
Tfull = tensor(Tfull);

[ lambda_est, Uest, maxeig, v ] = tensor_eig( Tfull ); % compute largest eigenvectors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [U,~,V] = svds(-gradf,1);
% % for smaller/non low rank matrices, use the next line instead of the above
% [U,~,V] = svd(-gradf,'econ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


atom1 = ktensor(1,v,v,v);
atom2=full(atom1);
atom = -atom2(:)/norm(atom2(:));


end

