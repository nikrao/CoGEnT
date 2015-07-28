function D = overDCTdict(n,L)
%


D = zeros(n,L);
D(:,1) = 1/sqrt(n);
for k = 2:L
  v = cos((0:n-1)*pi*(k-1)/L)';
  v = v-mean(v);
  D(:,k) = v/norm(v);
  if mod(k,100)==0
  fprintf('.')
  end
end