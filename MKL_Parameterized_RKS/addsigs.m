function [signew] = addsigs(sigold)


temp = sigold(2:end);
temp(end+1) = 0;

cs = (temp+sigold)/2;

signew = cs(1:end-1);

end