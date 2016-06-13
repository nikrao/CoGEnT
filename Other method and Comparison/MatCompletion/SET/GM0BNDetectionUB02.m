function [Uc,UcJ,h,Flag] = GM0BNDetectionUB02(X0,IndexM,UB,UBp,Uc)

% find aligned Uc and h
VS = Factorization01(X0,IndexM,UB*Uc,[]);
X_hat = ((UB*Uc)*VS').*IndexM;
gF = -2*UB'*(X0-X_hat)*VS;
gF = gF - Uc*Uc'*gF;
[UH SH VH] = svd(gF,0);
h = -UH(:,1);
Uc = Uc*VH;

% find angles
[tp,to,colv] = GM1dAnglesUB02(X0,IndexM,UB,UBp,Uc,h);

% detect 0BNs
tpv = tp(colv);
tov = to(colv);
cols_good = find(tpv < tov);
if isempty(cols_good) 
    UcJ = Uc;
    Flag = 0;
    return; 
end
% find the minimum tp such that there exist tos smaller than it. 
tom = repmat(tov(:),1,length(cols_good));
tpv = tpv(cols_good);
tpv = tpv(:);
tpm = repmat(tpv',length(colv),1);
i0BN = sum( tom < tpm,1 );
i0BN = i0BN>0;
if sum(i0BN) == 0
    UcJ = Uc;
    Flag = 0;
    return;
end
tpmin = min(tpv(i0BN));
tm = max( tov(tov<tpmin) );

% generate UcJ
Uct = Uc;
Uct(:,1) = Uc(:,1)*cos(tm)+h*sin(tm);
ht = -Uc(:,1)*sin(tm)+h*cos(tm);
VS = Factorization01(X0,IndexM,UB*Uct,[]);
X_hat = ((UB*Uct)*VS').*IndexM;
gFt = -2*UB'*(X0-X_hat)*VS;
corr = trace( [ht zeros(size(Uc,1),size(Uc,2)-1)]'*gFt*VH );
if corr == -1
    UcJ = Uct;
    Flag = 1;
else
    UcJ = Uc;
    Flag = 0;
end
