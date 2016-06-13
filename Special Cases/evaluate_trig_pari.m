function fval = evaluate_trig_pari(x,f,S, keepitreal)

n = length(x);

fval = zeros(length(f),1);
for j = 1:n
    
    if keepitreal
       fval = fval +  x(j)*cos(2*pi*f*S(j));
    else
        fval = fval + conj(x(j))*exp(1i*2*pi*f*S(j));
    end
end
fval = fval*sqrt(1/n);
fval = abs(fval);
 

end