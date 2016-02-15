function Immi = Immi(A, Cfull, varsizes)
% calculate redundancy between a set of Gaussian sources
% from minimum mutual information
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

chS = chol(Cs);
HS = sum(log(diag(chS))); % + 0.5*Nvary*log(2*pi*exp(1));
I = zeros(1,NA);
for ai=1:NA
    % use closed form expression
    chA = chol(AC(ai).Ca);
    chAS = chol(AC(ai).Cas);
    % normalisations cancel for information
    HA = sum(log(diag(chA))); % + 0.5*Nvarx*log(2*pi*exp(1));
    HAS = sum(log(diag(chAS))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));
    I(ai) = (HA + HS - HAS) / log(2);
end

if NA==1
    Immi = I(1);
end

if NA==2
    Immi = min(I);
end
