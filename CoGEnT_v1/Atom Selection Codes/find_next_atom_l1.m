function atom = find_next_atom_l1(gradf)

%find the next atom for l1 constrained problem
s = sign(-gradf);
[~,idx] = max(abs(gradf));
atom = sparse(zeros(size(gradf)));
atom(idx) = s(idx);

end