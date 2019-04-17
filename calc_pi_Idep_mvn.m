function Id = calc_pi_Idep_mvn(C, varsizes);
% calc_pi_Idep_mvn(C, varsizes)
% Calculate bivariate PID using Idep approach (James et al. 2017) for
% Gaussians (Kay & Ince 2018). 
% C is the full joint covariance matrix with variables in order
% X0 (first predictor), X1 (second predictor), S (target)
% varsizes is a length 3 vector containing the dimensionality of each of
% the above.
% Returns a length 4 vector containing
% [Red Unq_1 Unq_2 Syn]
%
% James et al. 2017 http://arxiv.org/abs/1709.06653
% Kay & Ince 2018 

vs = varsizes;
if length(vs)~=3
    error('calc_pi_Idep_mvn: only bivariate PID supported')
end
if sum(vs)~=size(C,1) || size(C,1)~=size(C,2)
    error('calc_pi_Idep_mvn: problem with variable specifications')
end

% variable indexes
xidx = 1:vs(1);
yidx = (vs(1)+1):(vs(1)+vs(2));
zidx = (vs(1)+vs(2)+1):(vs(1)+vs(2)+vs(3));

% extract blockwise covariance components
Cxx = C(xidx,xidx);
Cyy = C(yidx,yidx);
Czz = C(zidx,zidx);
Cxy = C(xidx,yidx);
Cxz = C(xidx,zidx);
Cyz = C(yidx,zidx);

% cholesky square root
chCxx = chol(Cxx);
chCyy = chol(Cyy);
chCzz = chol(Czz);

% Kay & Ince eq. D2
P = pinv(chCxx)'*Cxy*pinv(chCyy);
Q = pinv(chCxx)'*Cxz*pinv(chCzz);
R = pinv(chCyy)'*Cyz*pinv(chCzz);

% standard mutual informations
mivs = [vs(1)+vs(2)  vs(3)];
IXY = gauss_mi(C,mivs);
IY = gauss_mi(C([yidx zidx],[yidx zidx]), vs(2:3));
Ct = zeros(vs(1)+vs(3));
tzidx = (vs(1)+1):(vs(1)+vs(3));
Ct(xidx,xidx) = Cxx;
Ct(tzidx,tzidx) = Czz;
Ct(xidx,tzidx) = Cxz;
Ct(tzidx,xidx) = Cxz';
IX = gauss_mi(Ct, [vs(1) vs(3)]);

halflog2det = @(X) sum(log2(diag(chol(X))));
% Dependency lattice edges (Kay & Ince Table 9)
b = IX;

inum = halflog2det(eye(vs(2)) - R*Q'*Q*R');
iden = halflog2det(eye(vs(3))-Q'*Q) + halflog2det(eye(vs(3))-R'*R);
i = inum - iden - IY;

knum = halflog2det(eye(vs(2)) - P'*P);
Ct = zeros(sum(vs));
Ct(xidx,xidx) = eye(vs(1));
Ct(yidx,yidx) = eye(vs(2));
Ct(zidx,zidx) = eye(vs(3));
Ct(xidx,yidx) = P;
Ct(yidx,xidx) = P';
Ct(xidx,zidx) = Q;
Ct(zidx,xidx) = Q';
Ct(yidx,zidx) = R;
Ct(zidx,yidx) = R';
kden = halflog2det(Ct);
k = knum - kden - IY;

% fill out PID
Xunq = min([b i k]);
Id = zeros(1,4);
Id(2) = Xunq;
Id(1) = IX - Xunq;
Id(3) = IY - Id(1);
Id(4) = IXY - sum(Id(1:3));





