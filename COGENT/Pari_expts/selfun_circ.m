function atom = selfun_circ( gradf )
%  Function to do the atom selection step
%  Find the atom with maximum inner product with gradf
%  Note that the info about the atom (square wave, period, etc. is pasted
%  in here).
n=50;

gradf=reshape(-gradf,n,n);

I=eye(n,n);
P=circulant(I(2,:));
C=P;
C(4,8)=1;
C=C+C';

for i=1:n
    P2=P^(i-1);
    Ctemp=P2*C*P2';
    val(i)=sum(sum(Ctemp.*gradf));
    %Ctemp.*gradf
end
[t,idx]=max(val);
atom2=P^(idx-1)*C*(P^(idx-1))';
atom2=atom2.*gradf;
atom2=atom2/norm(atom2,'fro');
atom=atom2(:);

end