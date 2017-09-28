function pid = PID_dep_mvn(C, varsizes);
% 2 predictor PID with Idep for Gaussians

vs = varsizes;
if length(vs)~=3
    error('only bivaraite PID supported')
end
if sum(vs)~=size(C,1) || size(C,1)~=size(C,2)
    error('problem with variable specifications')
end

mivs = [vs(1)+vs(2)  vs(3)];

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

Cdiag = zeros(size(C));
Cdiag(xidx,xidx) = Cxx;
Cdiag(yidx,yidx) = Cyy;
Cdiag(zidx,zidx) = Czz;

IXY = gauss_mi(C,mivs);
IY = gauss_mi(C([yidx zidx]:[yidx zidx]), vs(2:3));
Ct = zeros(vs(1)+vs(3));
Ct(xidx,xidx) = Cxx;
Ct(zidx,zidx) = Czz;
Ct(xidx,zidx) = Cxz;
Ct(zidx,xidx) = Cxz';
IX = gauss_mi(Ct, [vs(1) vs(3)]);


% unique information in X
Xdelta = zeros(1,4);
Ctest_without = Cdiag;
Ctest_with = Ctest_without;
Ctest_with(xidx,zidx) = Cxz;
Ctest_with(zidx,xidx) = Cxz';
Xdelta(1) = gauss_mi(Ctest_with,mivs) - gauss_mi(Ctest_without,mivs);

Ctest_without = Cdiag;
Ctest_without(xidx,yidx) = Cxy;
Ctest_without(yidx,xidx) = Cxy';
Ctest_with = Ctest_without;
Ctest_with(xidx,zidx) = Cxz;
Ctest_with(zidx,xidx) = Cxz';
Xdelta(2) = gauss_mi(Ctest_with,mivs) - gauss_mi(Ctest_without,mivs);

Ctest_without = Cdiag;
Ctest_without(yidx,zidx) = Cyz;
Ctest_without(zidx,yidx) = Cyz';
Ctest_with = Ctest_without;
Ctest_with(xidx,zidx) = Cxz;
Ctest_with(zidx,xidx) = Cxz';
Xdelta(3) = gauss_mi(Ctest_with,mivs) - gauss_mi(Ctest_without,mivs);

Ctest_without = Cdiag;
Ctest_without(xidx,yidx) = Cxy;
Ctest_without(yidx,xidx) = Cxy';
Ctest_without(yidx,zidx) = Cyz;
Ctest_without(zidx,yidx) = Cyz';
Ctest_with = Ctest_without;
Ctest_with(xidx,zidx) = Cxz;
Ctest_with(zidx,xidx) = Cxz';
Xdelta(4) = gauss_mi(Ctest_with,mivs) - gauss_mi(Ctest_without,mivs);

Xunq = min(Xdelta);

% fill out PID
Id = zeros(1,4);
Id(2) = Xunq;
Id(1) = IX - Xunq;
Id(3) = IY - Id(1);
Id(4) = IXY - sum(Id(1:3));








