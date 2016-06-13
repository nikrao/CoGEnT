clear; close all;
randn('state',0); rand('state',0);
m = 9; n = 9; r = 2; 
SampRateV = [3:3:81]/(m*n);
srN = length(SampRateV);
tol = 1e-6;
n_o_realization = 6;


succRate = zeros(srN,n_o_realization);
for sn = 1:srN
    SampRate = SampRateV(sn);
    for nn = 1:n_o_realization
        [gU,gS,gV,b,Index,IndexRowCol]=LowRankMatrix03(m,n,r,SampRate);
        IndexM = repmat( 2==1,m,n );
        IndexM(Index) = repmat( 1==1,length(Index),1 );
        X0 = (gU*gS*gV').*IndexM;
        X0n2 = norm(X0,'fro')^2;

%         fprintf('\n sn=%d (SampRate=%f): nn = %d : \n',sn,SampRate,nn);
        [U,S,V,Xr_norm] = MatrixSET02(X0,IndexM,r,tol);
        if min(Xr_norm) < X0n2*tol
            succRate(sn,nn) = 1;
        end
        
        fprintf('.');
    end
    
    fprintf('\n');
end
disp(mean(succRate,2));
% save data_9_9_2 succRate;