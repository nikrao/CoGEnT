function [UB,UBp] = SpaceSplit03(X0,IndexM,U)
% split space: use best-matched columns

[m,n] = size(X0);

rLeft = 1:m;
UL = U;
UBd2 = 0;
UB = zeros(m,UBd2);
UBp = zeros(m,UBd2);
while 1
    % find columns that have length > r
    r = size(UL,2);
    ColNZN = sum(IndexM(rLeft,:),1);
    Colv = find( ColNZN>r );
    ColvLen = length( Colv );
    % if there is no such column existing
    if ColvLen == 0
        UB = [UB zeros(m,length(rLeft))];
        UB(rLeft,UBd2+1:UBd2+length(rLeft)) = eye( length(rLeft) );
        UBp = [UBp zeros(m,length(rLeft))];
        UBp(rLeft,UBd2+1:UBd2+length(rLeft)) = eye( length(rLeft) ); 
        UBd2 = UBd2+length(rLeft);
        break;
    end
    
    % find out how columns fit the UL
    X0L = X0(rLeft,Colv);
    X0Ln2 = sum( X0L.*X0L,1 );
    IndexML = IndexM(rLeft,Colv);
    VS = Factorization01(X0L,IndexML,UL,[]);
    XL_hat = UL*VS';
    XLr = zeros(size(X0L));
    XLr(IndexML) = X0L(IndexML)-XL_hat(IndexML);
    XLrn2 = sum( XLr.*XLr,1 );
    X0n2 = zeros(1,n);
    X0n2( Colv ) = X0Ln2;
    Xrn2 = zeros(1,n);
    Xrn2( Colv ) = XLrn2;
    0;
    
    % choose a column that has the best grade
    col_valid = find(X0n2>0);
    ratios = ones(1,n);
    ratios(col_valid) = Xrn2(col_valid)./X0n2(col_valid)./ColNZN(col_valid);
    [r_min,cn] = min(ratios);
    
    % form a block
    rnz_p= find(IndexM(rLeft,cn)==1);
    Uadd = orth(UL(rnz_p,:));
    Uadd_d2 = size(Uadd,2);
    rnz_p = rLeft(rnz_p);
    UB = [UB zeros(m,Uadd_d2)];
    UB( rnz_p,UBd2+1:UBd2+Uadd_d2 ) = Uadd;
    UBp = [UBp zeros(m,Uadd_d2)];
    UBp( rnz_p,UBd2+1:UBd2+Uadd_d2 ) = ones(length(rnz_p),Uadd_d2);
    UBd2 = UBd2+Uadd_d2;
    
    rLeft = setdiff(rLeft,rnz_p);
    UL = orth( U(rLeft,:) );
    
    if isempty(rLeft) 
        break;
    end
end
0;