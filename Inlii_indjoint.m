function Inli = Inlii(A, Pjoint)
% calculate redundancy from negative local interaction information
% A - cell array of elements
% Pjoint - full joint distribution

s = size(Pjoint);
Sm = s(end); % number of target values
Nx = length(s)-1; % number of dependent variables
vars = 1:Nx;

NA = length(A);
if NA>3
    error('Inli: only 3 elements supported')
end
Pele(NA).Pa = []; % intialize struct
Am = zeros(1,NA); % number of symbols in each element

% Ps
Ps = Pjoint;
for xi=1:Nx
    Ps = squeeze(sum(Ps,1));
end

% sort elements 
A = cellfun(@sort, A, 'Unif',false);
% build distributions for each element
for ai=1:NA
    thsA = A{ai};
    Nv = length(thsA);
    % vars to sum over
    sumover = setdiff(vars, thsA);
    Pas = Pjoint;
    for ii=1:length(sumover)
        Pas = sum(Pas, sumover(ii));
    end
    % joint distribution P(a,s)
    Pas = squeeze(Pas);
    % target first axis to collapse over non-target axes
%     Pas = permute(Pas, [Nv+1 1:Nv]);
%     Pas = Pas(:,:)';
    s = size(Pas);
    Pas = reshape(Pas, [prod(s(1:end-1)) s(end)]);
    Pele(ai).Pas = Pas;
    % unconditional distribution P(a)
    Pele(ai).Pa = squeeze(sum(Pas,2));
    Pele(ai).Pacs = bsxfun(@rdivide, Pele(ai).Pas, Ps);
    Am(ai) = size(Pas,1);
end

% build pairwise joint element distributions
if NA>1
    pairs = nchoosek(1:NA,2);
    Npair = size(pairs,1);
    Ppair(Npair).Paa = []; % intialize struct
    for pi=1:Npair
        p1 = pairs(pi,1);
        p2 = pairs(pi,2);
        Ppair(pi).Paacs = zeros(Am(p1),Am(p2),Sm);
        for si=1:Sm
            Paacs(:,:,si) = Pele(p1).Pacs(:,si) * Pele(p2).Pacs(:,si)';
        end
        Ppair(pi).Paacs = Paacs;
        Paas = bsxfun(@times, Paacs, reshape(Ps,[1 1 Sm]));
        Ppair(pi).Paas = Paas;
        Ppair(pi).Paa = nansum(Paas,3);
    end
end

% build triplewise joint element distributions
Paaas = cell(1,NA);
if NA==3
    Paaacs = zeros(Am(1),Am(2),Am(3),Sm);
    for si=1:Sm
        for a3=1:Am(3)
            for a2=1:Am(2)
                for a1=1:Am(1)
                    Paaacs(a1,a2,a3,si) = Pele(1).Pacs(a1,si) * Pele(2).Pacs(a2,si) * Pele(3).Pacs(a3,si);
                end
            end
        end
    end
    Ptrip(1).Paaacs = Paaacs;
    Paaas = bsxfun(@times, Paaacs, reshape(Ps,[1 1 1 Sm]));
    Ptrip(1).Paaas = Paaas;
    Ptrip(1).Paaa = nansum(Paaas,4);
end

% pointwise interaction information
pii = zeros([Am Sm]);
if NA==1
    for a1=1:Am(1)
        for si=1:Sm
            % local interaction information
            % = neg local mutual information
            num = Pele(1).Pa(a1) * Ps(si);
            den = Pele(1).Pas(a1,si);
            if den==0
                ii = 0;
            else
                ii = log2(num ./ den);
            end
            pii(a1,si) = ii;
        end
    end
    pii = Pele(1).Pas .* pii;
elseif NA==2
    for a1=1:Am(1)
        for a2=1:Am(2)
            for si=1:Sm
                num = Pele(1).Pa(a1) * Pele(2).Pa(a2) * Ps(si) * Ppair(1).Paas(a1,a2,si);
                den = Pele(1).Pas(a1,si) * Pele(2).Pas(a2,si) * Ppair(1).Paa(a1,a2);
                ii12 = log2(num ./ den);
%                 if num>0
%                     dsj = log2( Ppair(1).Paas(a1,a2,si) / (Ppair(1).Paa(a1,a2)*Ps(si)) );
%                     ds1 = log2( Pele(1).Pas(a1,si) ./ (Pele(1).Pa(a1)*Ps(si)) );
%                     ds2 = log2( Pele(2).Pas(a2,si) ./ (Pele(2).Pa(a2)*Ps(si)) );
%                     keyboard
%                 end
                
                num = Pele(1).Pa(a1) * Ps(si);
                den = Pele(1).Pas(a1,si);
                if den==0
                    ii1 = 0;
                else
                    ii1 = log2(num ./ den);
                end
                
                num = Pele(2).Pa(a2) * Ps(si);
                den = Pele(2).Pas(a2,si);
                if den==0
                    ii2 = 0;
                else
                    ii2 = log2(num ./ den);
                end
                
%                 pii(a1,a2,si) = nanmax([ii12 ii1 ii2]);
                pii(a1,a2,si) = ii12;
            end
        end
    end
    pii = Ppair(1).Paas .* pii;
elseif NA==3
    for a1=1:Am(1)
        for a2=1:Am(2)
            for a3=1:Am(3)
                for si=1:Sm
                    num = Pele(1).Pa(a1) * Pele(2).Pa(a2) * Pele(3).Pa(a3) * Ps(si);
                    num = num * Ppair(1).Paas(a1,a2,si) * Ppair(2).Paas(a1,a3,si) * Ppair(3).Paas(a2,a3,si);
                    num = num * Ptrip(1).Paaa(a1,a2,a3);

                    den = Pele(1).Pas(a1,si) * Pele(2).Pas(a2,si) * Pele(3).Pas(a3,si);
                    den = den * Ppair(1).Paa(a1,a2) * Ppair(2).Paa(a1,a3) * Ppair(3).Paa(a2,a3);
                    den = den * Ptrip(1).Paaas(a1,a2,a3,si);

                    ii123 = log2(num ./ den);

                    % pair(1) = 1 2
                    num = Pele(1).Pa(a1) * Pele(2).Pa(a2) * Ps(si) * Ppair(1).Paas(a1,a2,si);
                    den = Pele(1).Pas(a1,si) * Pele(2).Pas(a2,si) * Ppair(1).Paa(a1,a2);
                    ii12 = log2(num ./ den);

                    % pair(2) = 1 3
                    num = Pele(1).Pa(a1) * Pele(3).Pa(a3) * Ps(si) * Ppair(2).Paas(a1,a3,si);
                    den = Pele(1).Pas(a1,si) * Pele(3).Pas(a3,si) * Ppair(2).Paa(a1,a3);
                    ii13 = log2(num ./ den);

                    % pair(3) = 2 3
                    num = Pele(2).Pa(a2) * Pele(3).Pa(a3) * Ps(si) * Ppair(3).Paas(a2,a3,si);
                    den = Pele(2).Pas(a2,si) * Pele(3).Pas(a3,si) * Ppair(3).Paa(a2,a3);
                    ii23 = log2(num ./ den);

                    num = Pele(1).Pa(a1) * Ps(si);
                    den = Pele(1).Pas(a1,si);
                    if den==0
                        ii1 = 0;
                    else
                        ii1 = log2(num ./ den);
                    end

                    num = Pele(2).Pa(a2) * Ps(si);
                    den = Pele(2).Pas(a2,si);
                    if den==0
                        ii2 = 0;
                    else
                        ii2 = log2(num ./ den);
                    end

                    num = Pele(3).Pa(a3) * Ps(si);
                    den = Pele(3).Pas(a3,si);
                    if den==0
                        ii3 = 0;
                    else
                        ii3 = log2(num ./ den);
                    end

%                     pii(a1,a2,a3,si) = nanmax([ii123 ii12 ii13 ii23 ii1 ii2 ii3]);
                    % max over sub-pairs enforces monoticity
%                     pii(a1,a2,a3,si) = nanmax([ii123 ii12 ii13 ii23]);
                    % direct interaction information (not monotonic)
                    pii(a1,a2,a3,si) = ii123;
%                     if nansum(nanmax([ii123 ii12 ii13 ii23])) ~= nansum(nanmax(ii123)) && Ptrip(1).Paaas(a1,a2,a3,si)~=0
%                         keyboard
%                     end

                end
            end
        end
    end
    pii = Ptrip(1).Paaas .* pii;
end
% pii(~isfinite(pii))=0;
% pii
locred = -nansum(pii(pii<0));
Inli = locred;


function Pnew = copy_var(P, var, newpos)
% form joint distribution with variable var copied to axis position newpos
s = size(P);
varM = s(var);
% size of new array
news = [s(1:newpos-1) varM s(newpos:end)];
Pnew = zeros(news);

subP = cell(1,ndims(P));
[subP{:}] = ind2sub(size(P),1:numel(P));

subPnew = [subP(1:newpos-1) subP(var) subP(newpos:end)];
indPnew = sub2ind(size(Pnew), subPnew{:});
Pnew(indPnew) = P(:);



