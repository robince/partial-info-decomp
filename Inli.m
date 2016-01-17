function Inli = Inli(A, Pjoint)
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
    Am(ai) = size(Pas,1);
end

% build pairwise joint element distributions
if NA>1
    pairs = nchoosek(1:NA,2);
    Npair = size(pairs,1);
    Ppair(Npair).Paa = []; % intialize struct
    for pi=1:Npair
        thsA = [A{pairs(pi,1)} A{pairs(pi,2)}];
        Nv = length(thsA);
%         A1 = 1:length(A{pairs(pi,1)});
        Nv1 = length(A{pairs(pi,1)});
%         A2 = (length(A{pairs(pi,1)})+1):length(thsA);
        Nv2 = length(A{pairs(pi,2)});
        sumover = setdiff(vars, thsA);
        Paas = Pjoint;
        for ii=1:length(sumover)
            Paas = sum(Paas, sumover(ii));
        end
        % joint distribution over all variables
        % in both pairs of elements
        Paas = squeeze(Paas);
        % collapse A1
        s = size(Paas);
        Paas = reshape(Paas, [prod(s(1:Nv1)) s(Nv1+1:end)]);
        % collapse A2
        s = size(Paas);
        Paas = reshape(Paas, [s(1) prod(s(2:end-1)) s(end)]);
        Ppair(pi).Paas = Paas;
        Ppair(pi).Paa = squeeze(sum(Paas,3));
    end
end

% build triplewise joint element distributions
if NA==3
    Paaas = Pjoint;
    s = size(Paaas);
    % collapse A1
    Nv1 = length(A{1});
    Paaas = reshape(Paaas, [prod(s(1:Nv1)) s(Nv1+1:end)]);
    % collapse A2
    s = size(Paaas);
    Nv2 = length(A{2});
    Paas = reshape(Paas, [s(1) prod(s(2:Nv2+1)) s(Nv2+2:end)]);
    % collapse A3
    s = size(Paaas);
    Paas = reshape(Paas, [s(1:2) prod(s(3:end-1)) s(end)]);
    Ptrip(1).Paaas = Paaas;
    Ptrip(1).Paaa = squeeze(sum(Paaas,4));
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
            pii(a1,si) = log2(num ./ den);
        end
    end
    pii = Pele(1).Pas .* pii;
elseif NA==2
    for a1=1:Am(1)
        for a2=1:Am(2)
            for si=1:Sm
                num = Pele(1).Pa(a1) * Pele(2).Pa(a2) * Ps(si) * Ppair(1).Paas(a1,a2,si);
                den = Pele(1).Pas(a1,si) * Pele(2).Pas(a2,si) * Ppair(1).Paa(a1,a2);
                pii(a1,a2,si) = log2(num ./ den);
            end
        end
    end
    pii = Ppair(1).Paas .* pii;
elseif NA==3
    for a1=1:Am(1)
        for a2=1:Am(2)
            for a3=1:Am(3)
                for si=1:Sm
                    num = Pele(1).Pa(a1) * Pele(2).Pa(a2) * Ps(si) * Ppair(1).Paas(a1,a2,si);
                    den = Pele(1).Pas(a1,si) * Pele(2).Pas(a2,si) * Ppair(1).Paa(a1,a2);
                    pii(a1,a2,a3,si) = log2(num ./ den);
                end
            end
        end
    end
    pii = Ptrip(1).Paaas .* pii;
end

locred = -nansum(pii(pii<0));
Inli = locred;
