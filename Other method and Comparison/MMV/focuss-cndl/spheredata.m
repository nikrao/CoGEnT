function [out] = spheredata(in)
%  spheredata  - Subtract the mean and filter with the local symmetrical
%                zero-hpase whitening filter.  See Bell and Sejnowski
%                "The 'Independent Components' of Natural Scenes are Edge
%                 Filters"
%                Based on some code in:
%                ftp://ftp.cnl.salk.edu/pub/tony/sep96.public
%
%
% in        - Input data, each row is a sample
%
% Returns
% out       - Sphered (whitened) data
%
%
%  JFM   1/11/2001
%  Rev:  1/20/2001

[rows cols] = size(in);

inmean = mean(in')'; 
incov = cov(in');

out = in - inmean*ones(1,cols);                  

wz = 2*incov^(-1/2);
out = real(wz*out);             

return;