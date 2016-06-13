function w = Ax(x,A,groupINDS)

u = x'*groupINDS;
w = A*u';

end