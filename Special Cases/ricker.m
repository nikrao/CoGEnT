function y = ricker(t,s)

% returns the ricker wavelet evaluation of t with std. dev. = s
sigma = repmat(s^2,length(t),1);

y = (1 - (t.^2)./sigma);
y = y.*exp(-(t.^2)./(2*sigma));
y = y*2/(sqrt(3*s)*pi^(1/4));
end