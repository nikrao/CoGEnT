%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tensor complement algorithm implement 
% Version: 1.0
% Time: 03/01/2009
% Reference: "Tensor Completion for Estimating Missing Values 
% in Visual Data", ICCV 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y, errList, D] = LRTC(Y, alpha, beta, gamma, mark, maxIter, initial)
%%%%%%%%%%%%%%%%%%%%%%%%%%
% min(X, Y, M1, M2, M3,... Mn): (alpha1||X_(1)-M1||^2 + alpha2||X_(2)-T2||^2 + alpha3||X_(3)-T3||^2 + ...)/2 + 
%               beta1||Y_(1)-M1||^2 + beta2||Y_(2)-M2|| + beta3||Y_(3)-M3|| + ...+
%               gamma1||M1||_tr + gamma2||M2||_tr + gamma3||M3||_tr + ....
%         s.t.  Y_(1-ind) = T_(1-ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%
errList = zeros(maxIter, 1);
if isempty(initial)
    Y(mark) = 0; %mean(Y(logical(1-ind)));
else
    Y(mark) = initial(mark);
end

dim = size(Y);
M = cell(ndims(Y), 1);
for i = 1:ndims(Y)
    M{i} = Y;
end
X = Y;
Xsum = zeros(dim);
Ysum = zeros(dim);

Csum = alpha + beta;
for k = 1:maxIter
    k
    Xsum = Xsum * 0;
    Ysum = Ysum * 0;
    for i = 1:ndims(Y)
        site = circshift(dim, [1-i, 1-i]);
        Mpro = alpha(i)/Csum(i) * X + beta(i)/Csum(i) * Y;
        [mid, D(i)] = Pro2TraceNorm(reshape(shiftdim(Mpro,i-1), dim(i), []),...
            gamma(i)/Csum(i));
        M{i} = shiftdim(reshape(mid, site), -i+1+ndims(Y));
        Xsum = Xsum + alpha(i)*M{i};
        Ysum = Ysum + beta(i)*M{i};
    end
    X = Xsum / (sum(alpha) + 1e-10);
    Ysum = Ysum / (sum(beta) + 1e-10);
    errList(k) = norm(Y(mark)-Ysum(mark));
    Y(mark) = Ysum(mark);
end