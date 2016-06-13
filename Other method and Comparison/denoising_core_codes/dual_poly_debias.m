function [x_debias,fs,coeffs, pts,polyval]= dual_poly_debias( dual_poly_coeff, y)
% x = DUAL_POLY_DEBIAS( dual_poly_coeff, y, mu ) Debias by pole localization of dual polynomial
% Inputs
%  dual_poly_coeff: Coefficients of the dual polynomial
%  observed_samples: noisy mixture of sinusoids
%  mu: tradeoff parameter for atomic norm regularizer
n = length(y);
grid_size = 2^22; % size of the grid for finding poles.  Could possibly go up to 2^24?
amp = grid_size/sqrt(n); % adjust dual polynomial amplitude to peak at 1

polyval=amp*abs(ifft(conj(dual_poly_coeff(:)),grid_size)); % values of the dual polynomial on the grid
pts = linspace(0,1,grid_size);
% find all places where the dual polynomial is near 1
pole_regions = find(polyval > 1 - 1e-3);
if isempty(pole_regions)
    fs = [];
    coeffs = [];
    x_debias = 0;
    return;
end
% find connected segments where the polynomial is near 1
break_points = [1;find(diff(pole_regions)>1)+1;length(pole_regions)];

% in each connected segment, find the point closest to 1 and call this a
% pole
for k=1:length(break_points)-1,
    lookup = pole_regions(break_points(k):(break_points(k+1)-1));
    [~,slot] = max(abs(polyval(lookup)));
    %[~,slot] = min(abs(polyval(lookup)-1));
    pole_est(k) = 2*pi*lookup(slot)/grid_size;
end

fs = pole_est/2/pi;
% make a basis consisting of complex sinusoids at each pole
est_basis = exp( 1i*(0:(n-1))'*pole_est );
% fit the coefficients via least squares
coeffs = est_basis\y ;
% return the debiased signal
x_debias = est_basis*coeffs;

end
