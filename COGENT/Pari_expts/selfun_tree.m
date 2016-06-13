function atom = selfun_tree( gradf )
%  Function to do the atom selection step
%  Find the atom with maximum inner product with gradf
%  Note that the info about the atom (square wave, period, etc. is pasted
%  in here).
n=50;

gradf=reshape(gradf,n,n);
[T,C]=UndirectedMaximumSpanningTree(-gradf);
atom1=T/norm(T,'fro');
atom=atom1(:);


end