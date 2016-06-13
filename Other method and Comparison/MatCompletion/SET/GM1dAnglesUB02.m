function [tp,to,colv] = GM1dAnglesUB02(X0,IndexM,UB,UBp,Uc,hc)

n = size(X0,2);
r = size(Uc,2);
U = UB*Uc;
h = UB*hc;
u = U(:,1);
ur = repmat(u,1,n).*IndexM;
hr = repmat(h,1,n).*IndexM;

% calculate projection residues, including ur, hr and Xr;
if r>1
    U0 = U(:,2:r);
    VS = Factorization01(X0,IndexM,U0,[]);
    X_hat = (U0*VS').*IndexM;
    Xr = X0 - X_hat;
    
    urVS = Factorization01(ur,IndexM,U0,[]);
    ur = ur - (U0*urVS').*IndexM;
    hrVS = Factorization01(hr,IndexM,U0,[]);
    hr = hr - (U0*hrVS').*IndexM;
else
    Xr = X0;
end

% Given a column, identify part of Uc and h
overlap = UBp'*IndexM;
overlap = (overlap>0);
colv = find( sum(overlap,1)>r );

% calculate to
cor_uy = u'*Xr;
cor_hy = h'*Xr;
to = atan2( cor_uy.*sign(cor_uy) , -cor_hy.*sign(cor_uy) );

% calculate tp
tp = zeros(1,n);
colsN = length(colv);
pcoeff = zeros(2,n);
cn = 1;
while cn <= colsN
    col = colv(cn);
    if cond([ur(:,col) hr(:,col)]) > 1e10
        colv = colv([1:cn-1 cn+1:colsN]);
        colsN = colsN-1;
    else
        pcoeff(:,col) = [ur(:,col) hr(:,col)]\Xr(:,col);
        cn = cn+1;
    end
end
tp(colv) = atan2( pcoeff(2,colv).*sign(pcoeff(2,colv)) , ...
    pcoeff(1,colv).*sign(pcoeff(2,colv)) );
