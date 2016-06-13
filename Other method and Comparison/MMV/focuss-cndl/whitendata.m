function [X, whiten, dewhiten] = whitendata(X, rdim)
%  whitendata  - Subtract the mean and filter with the local symmetrical
%                zero-phase whitening filter.  
%
%                Based on some code at:
%                http://www.cis.hut.fi/projects/ica/imageica/
%
%
%
% in        - Input data, each column is a sample (wide matrix)
% whiten    - Whitening matrix 
% dewhiten  - Dewhitening matrix
% rdim      - Reducing to dimension rdim (if non-empty)
%
% Returns
% out       - Sphered (whitened) data
%
%
%  JFM   2/27/2001 (written by Patrik Hoyer, 9/1999)
%  Rev:  2/27/2001


%----------------------------------------------------------------------
% Subtract local mean gray-scale value from each patch
%----------------------------------------------------------------------

fprintf('Subtracting local mean...\n');
X = X-ones(size(X,1),1)*mean(X);

%----------------------------------------------------------------------
% Reduce the dimension and whiten at the same time!
%----------------------------------------------------------------------

% Calculate the eigenvalues and eigenvectors of covariance matrix.
fprintf ('Calculating covariance...\n');
covarianceMatrix = X*X'/size(X,2);
[E, D] = eig(covarianceMatrix);

% Sort the eigenvalues and select subset, and whiten

if(~isempty(rdim))
    fprintf('Reducing dimensionality (to %d) and whitening...\n', rdim);
    [dummy,order] = sort(diag(-D));
    E = E(:,order(1:rdim));
    d = diag(D); 
    d = real(d.^(-0.5));
    D = diag(d(order(1:rdim)));
end


%X = D*E'*X;

whiten = D*E';
dewhiten = E*D^(-1);

X = whiten * X;

return;

