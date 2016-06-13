function [Uc,VS,tv,fval] = GMLineMoveSplit1d01(X0,IndexM,U0,Uc,stepN)
% use golden ratio to move Uc
% calculate dynamic stepn

gRatio = (sqrt(5)-1)/2;
gRext = gRatio/(1-gRatio);
tvlen = ceil(9*log(10)/log(gRext))+stepN+1;

% calculate fval at U0
[m,n] = size(X0);
VS0 = Factorization01(X0,IndexM,U0*Uc,[]);
Xr = zeros(m,n);
X_hat = (U0*Uc)*VS0';
Xr(IndexM) = X0(IndexM) - X_hat(IndexM);
Xrn2 = sum( Xr(IndexM).*Xr(IndexM) );

% save data
fval = zeros(1,tvlen);
tv = zeros(1,tvlen);
tn = 1;
fval(tn) = Xrn2;
tv(tn) = 0;
tn = tn+1;

% calculate the gradient at U
dF = -2*U0'*Xr*VS0;
dF = dF - Uc*Uc'*dF;
[UdF SdF VdF] = svd(-dF,0);
SdF = diag(SdF);
% no need for update
if norm(SdF)^2/norm(X0(IndexM))^2 > 1e-10
    SdF = SdF * (pi/2/SdF(1));
    SdF(2:length(SdF)) = repmat(0,length(SdF)-1,1);
else
    VS = VS0;
    tv = tv(1:tn-1);
    fval = fval(1:tn-1);
    return;
end

% prepare for search, find the correct region of t, i.e., find the tmax
t = 1e-9;
while t <= 1
    Uct = [Uc*VdF UdF] * [diag(cos(SdF*t)); diag(sin(SdF*t))] * VdF';
    % fval at t
    VSt = Factorization01(X0,IndexM,U0*Uct,[]);
    X_hat = (U0*Uct)*VSt';
    Xr(IndexM) = X0(IndexM) - X_hat(IndexM);
    Xrn2 = sum(Xr(IndexM).*Xr(IndexM));
    % save data
    fval(tn) = Xrn2;
    tv(tn) = t;
    tn = tn+1;
    
    if fval(tn-1) < fval(tn-2)
        t = t*gRext;
    elseif t == 1e-9
        VS = VS0;
        tv = tv(1:tn-1);
        fval = fval(1:tn-1);
        return;
    else
        break;
    end
end
if t<=1
    tmax = t;
else
    tmax = tv(tn-1);
end

% prepare for search
fvalt = zeros(1,4); 
tvt(1)=tv(tn-3); tvt(2)=tv(tn-2); tvt(4)=tmax;
fvalt(1)=fval(tn-3); fvalt(2)=fval(tn-2); fvalt(4)=fval(tn-1);

stepn = 1;
tpos = 3;
tvt(3)=tvt(1)+gRatio*(tvt(4)-tvt(1));
while stepn <= stepN
    if tvt(4)-tvt(1) < 2e-10
        break;
    end
    
    t = tvt(tpos);
    Uct = [Uc*VdF UdF] * [diag(cos(SdF*t)); diag(sin(SdF*t))] * VdF';
    % fval at t
    VSt = Factorization01(X0,IndexM,U0*Uct,[]);
    X_hat = (U0*Uct)*VSt';
    Xr(IndexM) = X0(IndexM) - X_hat(IndexM);
    Xrn2 = sum(Xr(IndexM).*Xr(IndexM));
    fvalt(tpos) = Xrn2;    
    % save data
    fval(tn) = Xrn2;
    tv(tn) = t;
    tn = tn+1;
    
    % find the next step
    if fvalt(1)>fvalt(2) && fvalt(2)>=fvalt(3)
        % throw away the first point
        tvt(1:2) = tvt(2:3);
        fvalt(1:2) = fvalt(2:3);
        tvt(3) = (tvt(4)-tvt(1))*gRatio + tvt(1);
        tpos = 3;
    else
        % throw away the last point
        tvt(3:4) = tvt(2:3);
        fvalt(3:4) = fvalt(2:3);
        tvt(2) = (tvt(3)-tvt(1))*gRatio + tvt(1);
        tpos = 2;
    end
    
    stepn = stepn+1;
end

tv = tv(1:tn-1);
fval = fval(1:tn-1);
[fvalmin,tp] = min(fval);
Uct = [Uc*VdF UdF] * [diag(cos(SdF*tv(tp))); diag(sin(SdF*tv(tp)))] * VdF';
Uc = Uct;
VS = Factorization01(X0,IndexM,U0*Uc,[]);
