function [bits] = entropy(x, nbins)
%  entropy      - Estimates the entropy of vector sources by using 
%                 histogram approximation to PDF.  See Lewicki:1999
%                 and Lewicki:2000
%
%
%
% JFM    5/9/2001
% Rev:   6/20/2002

[rows, cols] = size(x);

nhist = hist(reshape(x, rows*cols, 1), nbins);
xdist = nhist ./ sum(nhist);  % sum(nhist) = N patterns * n coefficients (cols of A or rows of x)
bits = 0;
%bits2 = 0;

for i = 1:nbins
    if(xdist(i) ~= 0)
        bits = bits + (nhist(i) / cols) * log2(xdist(i));
       % bits2 = bits2 + xdist(i) * log2(xdist(i));
    end
end
bits = -bits;
%bits2 = -bits2;
 