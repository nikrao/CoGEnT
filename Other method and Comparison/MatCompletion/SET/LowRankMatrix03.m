function [U,S,V,b,Index,IndexRowCol]=LowRankMatrix03(m,n,r,SampRate)
% produce the low rank matrix X and the corresponding index matrix
% X = U*S*V';
U = randn(m,r); 
V = randn(n,r); 
S = randn(r,r);

% We want to form the matrix X = U*S*V'. We use some tricks to avoid
% golobal SVD
U0 = orth(U);
V0 = orth(V);
S0 = (U0'*U)*S*(V0'*V)';
[Ut S Vt] = svd(S0);
U = U0*Ut;
V = V0*Vt;
X = U*S*V';

% take samples
SampleN = ceil(SampRate*m*n);
SampleN = min(m*n,SampleN);
Index = randsample(m*n,SampleN);
b = X(Index);

IndexRowCol = zeros(SampleN,2);
IndexRowCol(:,1) = mod(Index-1,m)+1;
IndexRowCol(:,2) = ceil(Index/m);
0;