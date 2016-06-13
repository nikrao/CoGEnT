function [x_ast_debiased,Tu,fs_ast,cs_ast,pts,polyval,x_ast,tau] = ast_denoise(observed,noise_std)
%   AST_DENOISE   Run AST using ADMM by dual polynomial localization
%     [X_AST_DEBIASED,TU,FS_AST,CS_AST,PTS,X_AST,TAU] = AST_DENOISE(OBSERVED,NOISE_STD,VARARGIN)
% 
%   Given samples from a mixture of exponentials and an estimate of the
%   standard deviation of noise, this code runs the AST algorithm to find
%   the frequencies and amplitudes of the mixture of exponentials
%   
%   It also returns
%   
%   x_ast_debiased - debiased signal
%   Tu       - Low rank toeplitz matrix for frequency localization
%   x_ast    - estimated signal
%   fs_ast   - estimated frequencies (using dual polynomial)
%   cs_ast   - estimated amplitudes (debiased)
%   pts      - frequency values for dual polynomial
%   polyval  - dual polynomial values
%   tau      - Tradeoff parameter use
%   
%   Use varargin to optionally specify the total absolute tolerance
%   and relative tolerance to ADMM. The default values are 1e-4 and 1e-5.
%   Either specify both, or nothing.

n = length(observed);
tau = sqrt(log(n)+log(4*pi*log(n)))*noise_std; % Tradeoff

[x_ast,dual_poly_coeffs,Tu] = ast_admm(observed,tau);

[x_ast_debiased,fs_ast,cs_ast, pts,polyval] = dual_poly_debias(dual_poly_coeffs,observed);

% if isempty(fs_ast)
%     x_ast = x;
% end
end
