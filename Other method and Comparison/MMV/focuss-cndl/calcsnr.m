function [snrmean, snrstdev, snrk] = calcsnr(x, xest)
%  calcsnr      - Calculates the signal to noise ratio (dB) between x and xest.
%                 Signals must have been arranged/sorted properly.
%
% function [snr,snrk] = calcsnr(x, xest)
%               x - input signal (assumed that each signal is in one column)
%               xest - estimated signal 
% Returns:
%               snr - Average SNR (10log10) over all signals
%               snrk - SNR for each column
%
%
% JFM    6/13/2001
% Rev:   8/29/2003

[n, N] = size(x);

for k = 1:N
    xk = x(:,k);
    if(norm(xk) ~= 0) 
        xk = xk / norm(xk);
    end
    
    xestk = xest(:,k);
    if(norm(xestk) ~= 0) 
        xestk = xestk / norm(xestk);
    end
    
    resid = xk - xestk;
    
    % Account for possible sign change of learned vector
    snrk_minus = 10 * log10( norm(xk)^2 / norm(resid)^2 );
    
    resid = xk + xestk;
    snrk_plus = 10 * log10( norm(xk)^2 / norm(resid)^2 );
    
    snrk(k) = max([snrk_plus, snrk_minus]);
end

snrmean = mean(snrk);
snrstdev = std(snrk);

