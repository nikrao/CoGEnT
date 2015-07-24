function [u,v] = find_next_atom_nucnorm_matrix(gradf)

%find the next atom for nuclear norm constrained problem


[U,~,V,~,~] = lansvd(-gradf,1,'L'); % compute largest singular vectors

V = V';
u = U(:,1); % column vector
v = V(1,:); % row vector
u = u/norm(u);
v = v/norm(v);
clear U V;

end