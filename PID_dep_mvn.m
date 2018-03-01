function Id = PID_dep_mvn(C, varsizes);
% 2 predictor PID with Idep for Gaussians

vs = varsizes;
if length(vs)~=3
    error('only bivariate PID supported')
end
if sum(vs)~=size(C,1) || size(C,1)~=size(C,2)
    error('problem with variable specifications')
end

mivs = [vs(1)+vs(2)  vs(3)];
IXY = gauss_mi(C,mivs);

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

chCxx = chol(Cxx);
chCyy = chol(Cyy);
chCzz = chol(Czz);

% Kay & Ince eq. D2
P = pinv(chCxx)*Cxy*pinv(chCyy);
Q = pinv(chCxx)*Cxz*pinv(chCzz);
R = pinv(chCyy)*Cyz*pinv(chCzz);

IY = gauss_mi(C([yidx zidx],[yidx zidx]), vs(2:3));
Ct = zeros(vs(1)+vs(3));
tzidx = (vs(1)+1):(vs(1)+vs(3));
Ct(xidx,xidx) = Cxx;
Ct(tzidx,tzidx) = Czz;
Ct(xidx,tzidx) = Cxz;
Ct(tzidx,xidx) = Cxz';
IX = gauss_mi(Ct, [vs(1) vs(3)]);

% Dependency lattice edges (Kay & Ince Table 9)
% I(Y; Z)
b = IX;

halflog2det = @(X) sum(log2(diag(chol(X))));

inum = halflog2det(eye(vs(1)) - R*Q'*Q*R');
iden = halflog2det(eye(vs(2))-Q'*Q) + halflog2det(eye(vs(2))-R'*R);
i = inum - iden - IY;

knum = halflog2det(eye(vs(1)) - P'*P);
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





