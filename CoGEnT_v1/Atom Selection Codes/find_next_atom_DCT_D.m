function next_atom = find_next_atom_DCT_D(gradf,D)

% function to find most correlated DCT atom
gradf = -gradf;
n = length(gradf);
vals = gradf'*D(:,2:end);
[~,inds] = max(abs(vals));
next_atom = D(:,inds+1)*sign(vals(inds));
end
