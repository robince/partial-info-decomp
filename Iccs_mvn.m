function Iccs = Iccs_mvn(A, Cfull, varsizes)
% calculate redundancy between a set of Gaussian sources
% from pointwise common change in surprise
% using pairwise marginal maxent solution
%
% A is cell array of sources
% Cfull is full covariance of system
% varsizes specifies the number of variables in each X_i and S (S last)


if sum(varsizes) ~= size(Cfull,1)
    error('wrong number of variables specified')
end
if length(varsizes)~=3 && length(varsizes)~=4
    error('only 2 or 3 predictors supported')
end

NA = length(A);
NVs = varsizes(end);
Nx = length(varsizes)-1;
NVx = varsizes(1:end-1);
varstart = cumsum(varsizes)+1;
varstart = [1 varstart(1:end-1)];

uniquevars = unique([A{:}]);

sidx = varstart(end):(varstart(end)+NVs-1);
Cs = Cfull(sidx,sidx);

% build Cax for each source
AC = [];
for ai=1:NA
    thsA = A{ai};
    aidxfull = {};
    aidx = {};
    thsvstart = 1;
    for vi=1:length(thsA)
        aidxfull{vi} = varstart(thsA(vi)):(varstart(thsA(vi))+NVx(thsA(vi))-1);
        thsL = length(aidxfull{vi});
        aidx{vi} = thsvstart:(thsvstart+thsL-1);
        thsvstart = thsvstart+thsL;
    end
    thsNv = length(cell2mat(aidx));
    Cas = zeros(thsNv+NVs);
    Ca = zeros(thsNv);
    
    % fill in blocks
    % diagonal
    for vi=1:length(thsA)
        Ca(aidx{vi},aidx{vi}) = Cfull(aidxfull{vi},aidxfull{vi}); 
    end
    % off diagonal
    for vi=1:length(thsA)
        for vj=1:length(thsA)
            if vi==vj
                continue
            end
            Ca(aidx{vi},aidx{vj}) = Cfull(aidxfull{vi},aidxfull{vj});
        end
    end
    
    Cas(1:thsNv,1:thsNv) = Ca;
    % joint with S
    % diagonal
    thssidx = thsNv+1:thsNv+NVs;
    Cas(thssidx,thssidx) = Cs;
    % off diagonal
    for vi=1:length(thsA)
        Cas(aidx{vi},thssidx) = Cfull(aidxfull{vi},sidx);
        Cas(thssidx,aidx{vi}) = Cfull(sidx,aidxfull{vi});
    end
    
    Casoff = Cas(1:thsNv,thssidx);
    CXYY1 = Casoff * pinv(Cs);
    Cacs = Ca - CXYY1*Cas(thssidx,1:thsNv);
    MacsF = CXYY1;
    
    AC(ai).Ca = Ca;
    AC(ai).Cas = Cas;
    AC(ai).Cacs = Cacs;
    AC(ai).Casoff = Casoff;
    AC(ai).MacsF = CXYY1;
    AC(ai).Nv = thsNv + NVs;
    AC(ai).ANv = thsNv;
    AC(ai).aidxfull = aidxfull;
    AC(ai).aidx = aidx;
end

% build Caas
if NA==2
    % NEED TO AVOID DUPLICATED VARS 
%     a12idxfull = cell2mat(cat(2,AC(1).aidxfull,AC(2).aidxfull));
%     thsNv = AC(1).ANv + AC(2).ANv + NVs;
%     ANv = AC(1).ANv + AC(2).ANv;
    
    thsA = unique(cat(2,A{1},A{2}));
    a12idxfull = {};
    aidx = {};
    thsvstart = 1;
    for vi=1:length(thsA)
        a12idxfull{vi} = varstart(thsA(vi)):(varstart(thsA(vi))+NVx(thsA(vi))-1);
        thsL = length(a12idxfull{vi});
        aidx{vi} = thsvstart:(thsvstart+thsL-1);
        thsvstart = thsvstart+thsL;
    end
    thsNv = length(cell2mat(aidx));
    a12idxfull = cell2mat(a12idxfull);
    % P2 == full gaussian covariance
    Caa = Cfull(a12idxfull,a12idxfull);
    Caas = Cfull([a12idxfull sidx],[a12idxfull sidx]);
    Caasoff = Cfull(a12idxfull,sidx);
    CXYY1 = Caasoff * pinv(Cs);
    Caacs = Caa - CXYY1*Cfull(sidx,a12idxfull);
    MaacsF = CXYY1;
    
    AAC(1).Caa = Caa;
    AAC(1).Caas = Caas;
    AAC(1).Caacs = Caacs;
    AAC(1).MaacsF = CXYY1;
    AAC(1).Nv = thsNv;  
    AAC(1).ANv = length(thsA);
    AAC(1).thsA = thsA;
    AAC(1).aidx = aidx;
    a1idx = [];
    for ai=1:length(A{1})
        a1idx = [a1idx aidx{thsA==ai}];
    end
    a2idx = [];
    for ai=1:length(A{2})
        a2idx = [a2idx aidx{thsA==ai}];
    end
    AAC(1).a1idx = a1idx;
    AAC(1).a2idx = a2idx;
        
end

% build Caas
if NA==3
    % extract required sources to form full three way A cov matrix
    
    thsA = unique(cat(2,A{1},A{2},A{3}));
    a123idxfull = {};
    aidx = {};
    thsvstart = 1;
    for vi=1:length(thsA)
        a123idxfull{vi} = varstart(thsA(vi)):(varstart(thsA(vi))+NVx(thsA(vi))-1);
        thsL = length(a123idxfull{vi});
        aidx{vi} = thsvstart:(thsvstart+thsL-1);
        thsvstart = thsvstart+thsL;
    end
    thsNv = length(cell2mat(aidx));
    a123idxfull = cell2mat(a123idxfull);
    
    % P2 == full gaussian covariance
    Caaa = Cfull(a123idxfull,a123idxfull);
    Caaas = Cfull([a123idxfull sidx],[a123idxfull sidx]);
    Caaasoff = Cfull(a123idxfull,sidx);
    CXYYY1 = Caaasoff * pinv(Cs);
    Caaacs = Caaa - CXYYY1*Cfull(sidx,a123idxfull);
    MaaacsF = CXYYY1;
        
    AAAC(1).Caaa = Caaa;
    AAAC(1).Caaas = Caaas;
    AAAC(1).Caaacs = Caaacs;
    AAAC(1).MaaacsF = CXYYY1;
    AAAC(1).Nv = thsNv;  
    AAAC(1).ANv = length(thsA);
    AAAC(1).thsA = thsA;
    AAAC(1).aidx = aidx;
    a1idx = [];
    for ai=1:length(A{1})
        a1idx = [a1idx aidx{thsA==ai}];
    end
    a2idx = [];
    for ai=1:length(A{2})
        a2idx = [a2idx aidx{thsA==ai}];
    end
    a3idx = [];
    for ai=1:length(A{3})
        a3idx = [a3idx aidx{thsA==ai}];
    end
    AAAC(1).a1idx = a1idx;
    AAAC(1).a2idx = a2idx;
    AAAC(1).a3idx = a2idx;
    a12idx = [];
    A12 = unique(cat(2,A{1},A{2}));
    for ai=1:length(A12)
        a12idx = [a12idx aidx{A12==ai}];
    end
    a13idx = [];
    A13 = unique(cat(2,A{1},A{3}));
    for ai=1:length(A13)
        a13idx = [a13idx aidx{A13==ai}];
    end
    a23idx = [];
    A23 = unique(cat(2,A{2},A{3}));
    for ai=1:length(A23)
        a23idx = [a23idx aidx{A23==ai}];
    end
    AAAC(1).a12idx = a12idx;
    AAAC(1).a13idx = a13idx;
    AAAC(1).a23idx = a23idx;

    % extract for each pair
    pairs = nchoosek(1:3,2);
    Npair = size(pairs,1);
    Ppair(Npair).Paa = []; % intialize struct
    for pi=1:Npair
        keepax = [pairs(pi,1) pairs(pi,2)];
        p1 = pairs(pi,1);
        p2 = pairs(pi,2);
        
        thsA = unique(cat(2,A{p1},A{p2}));
        a12idxfull = {};
        aidx = {};
        thsvstart = 1;
        for vi=1:length(thsA)
            a12idxfull{vi} = varstart(thsA(vi)):(varstart(thsA(vi))+NVx(thsA(vi))-1);
            thsL = length(a12idxfull{vi});
            aidx{vi} = thsvstart:(thsvstart+thsL-1);
            thsvstart = thsvstart+thsL;
        end
        thsNv = length(cell2mat(aidx));
        a12idxfull = cell2mat(a12idxfull);
        
        Caa = Cfull(a12idxfull,a12idxfull);
        Caas = Cfull([a12idxfull sidx],[a12idxfull sidx]);
        Caasoff = Cfull(a12idxfull,sidx);
        CXYY1 = Caasoff * pinv(Cs);
        Caacs = Caa - CXYY1*Cfull(sidx,a12idxfull);
        MaacsF = CXYY1;
        
        AAC(pi).Caa = Caa;
        AAC(pi).Caas = Caas;
        AAC(pi).Caacs = Caacs;
        AAC(pi).MaacsF = CXYY1;
        AAC(pi).Nv = thsNv;
        AAC(pi).ANv = length(thsA);
    end
end

%
% Calculate Iccs 
%

% Monte Carlo Integration for Iccs
% 100,000 works well for 5d space.
% 10,000 might be ok but some variance
Nmc = 1000000;

if NA==1
    % use closed form expression
    chA = chol(AC(1).Ca);
    chS = chol(Cs);
    chAS = chol(AC(1).Cas);
    % normalisations cancel for information
    HA = sum(log(diag(chA))); % + 0.5*Nvarx*log(2*pi*exp(1));
    HS = sum(log(diag(chS))); % + 0.5*Nvary*log(2*pi*exp(1));
    HAS = sum(log(diag(chAS))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));
    Iccs = (HA + HS - HAS);
end

if NA==2
    a1idx = AAC(1).a1idx;
    a2idx = AAC(1).a2idx;
    a12idx = 1:AAC(1).ANv;
    thssidx = AAC(1).ANv+1:AAC(1).ANv+NVs;

    % sample from multivariate gaussian we want expectation over.
    mcx = mvnrnd(zeros(1,AAC(1).Nv+NVs), AAC(1).Caas, Nmc);
    
    px = logmvnpdf(mcx(:,a1idx),zeros(1,AC(1).ANv),AC(1).Ca);
    pxcs = logmvnpdf(mcx(:,a1idx), (AC(1).MacsF*mcx(:,thssidx)')', AC(1).Cacs);
    
    py = logmvnpdf(mcx(:,a2idx),zeros(1,AC(2).ANv),AC(2).Ca);
    pycs = logmvnpdf(mcx(:,a2idx), (AC(2).MacsF*mcx(:,thssidx)')', AC(2).Cacs);
    
    pxy = logmvnpdf(mcx(:,a12idx), zeros(1,AAC(1).ANv), AAC(1).Caa);
    pxycs = logmvnpdf(mcx(:,a12idx), (AAC(1).MaacsF*mcx(:,thssidx)')', AAC(1).Caacs);
    
    dhx = pxcs - px;
    dhy = pycs - py;
    dhxy = pxycs - pxy;
    
    coi = dhx + dhy - dhxy;
    keep = sign(dhx) == sign(coi) & sign(dhx) == sign(dhy) & sign(dhx) == sign(dhxy);
    coi(~keep) = 0;
    Iccs = nanmean(coi);
end

if NA==3
    a1idx = AAAC(1).a1idx;
    a2idx = AAAC(1).a2idx;
    a3idx = AAAC(1).a3idx;
    a12idx = AAAC(1).a12idx;
    a13idx = AAAC(1).a13idx;
    a23idx = AAAC(1).a23idx;
    a123idx = 1:AAAC(1).ANv;
    thssidx = AAAC(1).ANv+1:AAAC(1).ANv+NVs;
    
    % sample from multivariate gaussian we want expectation over.
    mcx = mvnrnd(zeros(1,AAAC(1).ANv+NVs), AAAC(1).Caaas, Nmc);
    
    p1 = logmvnpdf(mcx(:,a1idx),zeros(1,AC(1).ANv),AC(1).Ca);
    p1cs = logmvnpdf(mcx(:,a1idx), (AC(1).MacsF*mcx(:,thssidx)')', AC(1).Cacs);
    
    p2 = logmvnpdf(mcx(:,a2idx),zeros(1,AC(2).ANv),AC(2).Ca);
    p2cs = logmvnpdf(mcx(:,a2idx), (AC(2).MacsF*mcx(:,thssidx)')', AC(2).Cacs);
    
    p3 = logmvnpdf(mcx(:,a3idx),zeros(1,AC(3).ANv),AC(3).Ca);
    p3cs = logmvnpdf(mcx(:,a3idx), (AC(3).MacsF*mcx(:,thssidx)')', AC(3).Cacs);
    
    % NB: hardcoding order output of nchoosek here, OK for 2 from 3 I
    % think. 
    p12 = logmvnpdf(mcx(:,a12idx), zeros(1,AAC(1).ANv), AAC(1).Caa);
    p12cs = logmvnpdf(mcx(:,a12idx), (AAC(1).MaacsF*mcx(:,thssidx)')', AAC(1).Caacs);
    
    p13 = logmvnpdf(mcx(:,a13idx), zeros(1,AAC(2).ANv), AAC(2).Caa);
    p13cs = logmvnpdf(mcx(:,a13idx), (AAC(2).MaacsF*mcx(:,thssidx)')', AAC(2).Caacs);
    
    p23 = logmvnpdf(mcx(:,a23idx), zeros(1,AAC(3).ANv), AAC(3).Caa);
    p23cs = logmvnpdf(mcx(:,a12idx), (AAC(1).MaacsF*mcx(:,thssidx)')', AAC(1).Caacs);
    
    
    p123 = logmvnpdf(mcx(:,a123idx), zeros(1,AAAC(1).ANv), AAAC(1).Caaa);
    p123cs = logmvnpdf(mcx(:,a123idx), (AAAC(1).MaaacsF*mcx(:,thssidx)')', AAAC(1).Caaacs);

    dh1 = p1cs - p1;
    dh2 = p2cs - p2;
    dh3 = p3cs - p3;
    dh12 = p12cs - p12;
    dh13 = p13cs - p13;
    dh23 = p23cs - p23;
    dh123 = p123cs - p123;
    
%     coi = dhx + dhy - dhxy;
    coi = dh1 + dh2 + dh3 - dh12 - dh13 - dh23 + dh123;
    keep = sign(dh1) == sign(coi) ...
        & sign(dh1) == sign(dh2) ...
        & sign(dh1) == sign(dh3) ...
        & sign(dh1) == sign(dh123);
    coi(~keep) = 0;
    Iccs = nanmean(coi);
end

% convert to bits
Iccs = Iccs ./ log(2);
    
