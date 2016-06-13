function [x_debiased,X,fs,cs,x] = cadzow_denoise(y, r)
% Alternately project onto Rank-r Grasmmanian manifold
% and space of Toeplitz matrices.
x = y;
max_iterations = 5000;
tol = 1e-10;
ratio = tol;

n = length(x);
X = T(x);

for iter=1:max_iterations
  if (ratio<tol), break; end
  % Alternate Projection using SVD, toeplitz_approx
  [U,D,V] = svd(X);
  V((r+1):end,(r+1):end)=0;
	d = diag(D);
  X = toeplitz_approx(U(:,1:r)*diag(d(1:r))*V(:,1:r)');
  % Update convergence parameter
  ratio = D(r+1,r+1)/D(r,r);
end

x = invT(x,n);
[x_debiased,fs,cs] = matrix_pencil_debias(invT(X,n), y, r);

end
