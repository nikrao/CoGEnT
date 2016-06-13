function atom = selfun_square( gradf )
%  Function to do the atom selection step
%  Find the atom with maximum inner product with gradf
%  Note that the info about the atom (square wave, period, etc. is pasted
%  in here).

N_grid=100;
N=100;
disc=linspace(0,2*pi,N_grid);
disc_freq=linspace(0,1,N_grid);
t=linspace(0,10,N);

for i=1:N_grid;
    for j=1:N_grid
    candidate_atom=sign(sin(disc_freq(j)*t-disc(i)))';
    inner_product(i,j)=gradf'*candidate_atom;
    end
end


[val, location]=max(inner_product(:));
[R,C] = ind2sub(size(inner_product),location);
atom=-sign(sin(disc_freq(C)*t-disc(R)))';

end