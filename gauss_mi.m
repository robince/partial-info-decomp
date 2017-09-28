function I = gauss_mi(C, varsizes)

% C = C + 0.1*eye(size(C));
vs = varsizes;
if length(vs)>2
    error('bivariate MI only')
end
if sum(vs) ~= size(C,1)
    error('variables incorrectly specified')
end

CX = C(1:vs(1),1:vs(1));
CY = C((vs(1)+1):end, (vs(1)+1):end);

% use closed form expression
chX = chol(CX);
chY = chol(CY);
chXY = chol(C);
% normalisations cancel for information
HX = sum(log(diag(chX))); % + 0.5*Nvarx*log(2*pi*exp(1));
HY = sum(log(diag(chY))); % + 0.5*Nvary*log(2*pi*exp(1));
HXY = sum(log(diag(chXY))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));
I = (HX + HY - HXY) / log(2);
