function Iccs = Iccs_mvn_P2(A, Cfull, varsizes)
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
if length(varsizes)~=3
    error('only 2 variables supported')
end

NA = length(A);
NVs = varsizes(end);
Nx = length(varsizes-1);
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
    AC(ai).Nv = thsNv;
end


if NA==1
    % use closed form expression
    chA = chol(AC(1).Ca);
    chS = chol(Cs);
    chAS = chol(AC(1).Cas);
    % normalisations cancel for information
    HA = sum(log(diag(chA))); % + 0.5*Nvarx*log(2*pi*exp(1));
    HS = sum(log(diag(chS))); % + 0.5*Nvary*log(2*pi*exp(1));
    HAS = sum(log(diag(chAS))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));
    Iccs = (HA + HS - HAS) / log(2);
end

if NA==2
    % Covariance for Pind(A1,A2)
    thsNv = AC(1).Nv + AC(2).Nv + NVs;
    ANv = AC(1).Nv + AC(2).Nv;
    a1idx = 1:AC(1).Nv;
    a2idx = AC(1).Nv+1:AC(1).Nv+AC(2).Nv;
    a12idx = 1:AC(1).Nv + AC(2).Nv;
    
    % P2 == full gaussian covariance
    C = Cfull(a12idx,a12idx);
    
    Caasoff = Cfull(a12idx,sidx);
    CXYY1 = Caasoff * pinv(Cs);
    Caacs = C - CXYY1*Cfull(sidx,a12idx);
    MaacsF = CXYY1;
    
    thssidx = AC(1).Nv+AC(2).Nv+1:AC(1).Nv+AC(2).Nv+NVs;

    % Monte Carlo Integration for Iccs
    % 100,000 works well for 5d space. 
    % 10,000 might be ok but some variance
    Nmc = 100000;
    
    % integrate over full space
    intD = AC(1).Nv + AC(2).Nv + NVs;

    % sample from multivariate gaussian we want expectation over.
    mcx = mvnrnd(zeros(1,ANv+NVs), Cfull, Nmc);
    
    px = logmvnpdf(mcx(:,a1idx),zeros(1,AC(1).Nv),AC(1).Ca);
    pxcs = logmvnpdf(mcx(:,a1idx), (AC(1).MacsF*mcx(:,thssidx)')', AC(1).Cacs);
    
    py = logmvnpdf(mcx(:,a2idx),zeros(1,AC(2).Nv),AC(2).Ca);
    pycs = logmvnpdf(mcx(:,a2idx), (AC(2).MacsF*mcx(:,thssidx)')', AC(2).Cacs);
    
    pxy = logmvnpdf(mcx(:,a12idx), zeros(1,ANv), C);
    pxycs = logmvnpdf(mcx(:,a12idx), (MaacsF*mcx(:,thssidx)')', Caacs);
    
    dhx = pxcs - px;
    dhy = pycs - py;
    dhxy = pxycs - pxy;
    
    lnii = dhx + dhy - dhxy;
    keep = sign(dhx) == sign(lnii) & sign(dhx) == sign(dhy);
    lnii(~keep) = 0;
   
    Iccs = nanmean(lnii);
    % convert to bits
    Iccs = Iccs ./ log(2);
end
