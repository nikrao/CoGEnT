function atom = selfun_sawtooth( gradf )
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
    candidate_atom=sawtooth(disc_freq(j)*t-disc(i))'/norm(sawtooth(disc_freq(j)*t-disc(i)));
    inner_product(i,j)=gradf'*candidate_atom;
    end
end


[val, location]=min(inner_product(:));
[R,C] = ind2sub(size(inner_product),location);
atom=sawtooth(disc_freq(C)*t-disc(R))'/norm(sawtooth(disc_freq(C)*t-disc(R)));

end