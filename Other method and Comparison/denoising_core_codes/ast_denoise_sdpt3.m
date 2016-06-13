function [ x_ast_debiased, Tu, fs_ast, cs_ast, pts, polyval, x_ast, tau ] = ast_denoise( observed, noise_std )
%AST_DENOISE Run AST using ADMM, debias by dual polynomial localization
%   Needs estimate of noise standard deviations
% Returns:
%   x_ast   - estimated signal
%   f_ast   - estimated frequencies (using dual polynomial)
%   c_ast   - estimated amplitudes (debiased)
%   pts     - frequency values for dual polynomial
%   polyval - dual polynomial values
n = length(observed);
tau = sqrt(log(n)+log(4*pi*log(n)))*noise_std; % Tradeoff
% tau = tau*(1+1/log(n));
% tau = sqrt(log(n))*noise_std;
[x_ast,dual_poly_coeffs,Tu] = ast_sdpt3(observed,tau);

[x_ast_debiased,fs_ast,cs_ast, pts,polyval] = dual_poly_debias(dual_poly_coeffs,observed);
% if isempty(fs_ast)
%     x_ast = x;
% end
end

