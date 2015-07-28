function atom = find_next_atom_linf(gradf)

%find the next atom for l-inf constrained problem

gradf = -gradf; % maximally correlate with this
atom = sign(gradf);
idx = find(atom == 0);
vec = 2*binornd(1,0.5,length(idx),1)-1;
atom(idx) = vec;

end