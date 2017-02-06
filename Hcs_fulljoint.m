function Hcs = Hcs_fulljoint(A, Pjoint)
% calculate redundant entropy between a set of sources
% from pointwise common surprisal
% using full joint 3-way distribution (no maxent)
%
% A - cell array of sources
% Pjoint - full joint distribution (S last axis)

isclosefp = @(a,b) abs(a - b) <= eps(max(abs(a), abs(b)));
s = size(Pjoint);
Nx = length(s); % number of dependent variables
vars = 1:Nx;

NA = length(A);
if NA>3
    error('Hcs: only 3 sources supported')
end
PA(NA).Pa = []; % intialize struct
Am = zeros(1,NA); % number of symbols in each element

% sort variables within each source
A = cellfun(@sort, A, 'Unif',false);

% build distributions for each source
for ai=1:NA
    thsA = A{ai};
    Nv = length(thsA);
    % vars to sum over
    sumover = setdiff(vars, thsA);
    Pa = Pjoint;
    for ii=1:length(sumover)
        Pa = sum(Pa, sumover(ii));
    end
    % distribution P(a)
    Pa = squeeze(Pa);
    Pa = Pa(:);
    PA(ai).Pa = Pa;
    Am(ai) = size(Pa,1);
end

% build pairwise joint element distributions
if NA>1
    pairs = nchoosek(1:NA,2);
    Npair = size(pairs,1);
    Ppair(Npair).Paa = []; % intialize struct
    for pi=1:Npair
        thsA = [A{pairs(pi,1)} A{pairs(pi,2)}];
        Nv = length(thsA);
        Nv1 = length(A{pairs(pi,1)});
        Nv2 = length(A{pairs(pi,2)});

        % collapse variables we don't need
        sumover = setdiff(vars, thsA);
        Paa = Pjoint;
        for ii=1:length(sumover)
            Paa = sum(Paa, sumover(ii));
        end
        Paa = squeeze(Paa);

        % reorder axes to match order of unique variables in this pair of
        % elements
        % order we want
        Aunq = unique(thsA,'stable');
        % order we have
        [Aunqsrt, Aunqsrtidx] = sort(Aunq);
        % invert order
        [~, Aidx] = sort(Aunqsrtidx);
%         Paa = permute(Paa, [Aidx length(Aunq)+1]);
        Paa = permute(Paa, [Aidx]);
        thsA = changem(thsA, 1:length(Aunq), Aunq);
        Aunq = unique(thsA, 'stable');

        % copy duplicate variables as required
        uniquevar_i = 1;
        for allvar_i=1:Nv
            if (uniquevar_i>length(Aunq)) || (thsA(allvar_i) ~= Aunq(uniquevar_i))
                % need to insert a duplicate variable
                var_needed = thsA(allvar_i);
                copy_from = find(thsA==var_needed,1);
                Paa = copy_var(Paa, copy_from, allvar_i);
            else
                % axis order is correct
                uniquevar_i = uniquevar_i + 1;
            end
        end

        % joint distribution over all variables
        % in both pairs of elements
        % now should have correct variable axis in correct order
        % collapse A1
        s = size(Paa);
        Paa = reshape(Paa, [prod(s(1:Nv1)) s(Nv1+1:end)]);
        % collapse A2
        s = size(Paa);
        Paa = reshape(Paa, [s(1) prod(s(2:end))]);
        Ppair(pi).Paa = Paa;
    end
end

% build triplewise joint element distributions
Paaa = cell(1,NA);
if NA==3
    thsA = [A{1} A{2} A{3}];
    Nv = length(thsA);
    Nv1 = length(A{1});
    Nv2 = length(A{2});
    Nv3 = length(A{3});

    % collapse variables we don't need
    sumover = setdiff(vars, thsA);
    Paaa = Pjoint;
    for ii=1:length(sumover)
        Paaa = sum(Paaa, sumover(ii));
    end
    Paaa = squeeze(Paaa);

    % reorder axes to match order of unique variables in this pair of
    % elements
    % order we want
    Aunq = unique(thsA,'stable');
    % order we have
    [Aunqsrt, Aunqsrtidx] = sort(Aunq);
    % invert order
    [~, Aidx] = sort(Aunqsrtidx);

    Paaas= permute(Paaa, [Aidx length(Aunq)+1]);
    thsA = changem(thsA, 1:length(Aunq), Aunq);
    Aunq = unique(thsA, 'stable');

    % copy duplicate variables as required
    uniquevar_i = 1;
    for allvar_i=1:Nv
        if (uniquevar_i>length(Aunq)) || (thsA(allvar_i) ~= Aunq(uniquevar_i))
            % need to insert a duplicate variable
            var_needed = thsA(allvar_i);
            copy_from = find(thsA==var_needed,1);
            Paaa = copy_var(Paaa, copy_from, allvar_i);
        else
            % axis order is correct
            uniquevar_i = uniquevar_i + 1;
        end
    end

    % joint distribution over all variables
    % now should have correct variable axes in correct order
    % collapse A1
    s = size(Paaa);
    Nv1 = length(A{1});
    Paaa = reshape(Paaa, [prod(s(1:Nv1)) s(Nv1+1:end)]);
    % collapse A2
    s = size(Paaa);
    Nv2 = length(A{2});
    Paaa = reshape(Paaa, [s(1) prod(s(2:Nv2+1)) s(Nv2+2:end)]);
    % collapse A3
    s = size(Paaa);
    Paaa = reshape(Paaa, [s(1:2) prod(s(3:end))]);
    Ptrip(1).Paaa = Paaa;
end

% pointwise co-information
cs = zeros([Am 1]);
if NA==1
    for a1=1:Am(1)
        s1 = -log2( PA(1).Pa(a1) );
        cs(a1) = s1;
    end
%     keyboard
    cs = PA(1).Pa .* cs;
elseif NA==2
    for a1=1:Am(1)
        for a2=1:Am(2)
            
            s1 = -log2( PA(1).Pa(a1) );
            s2 = -log2( PA(2).Pa(a2) );
            sj = -log2( Ppair(1).Paa(a1,a2) );
            
            % local co-information (entropy overlap)
            i = (s1 + s2 - sj);
            % local entropy always positive
            % if local co-information positive then have overlap
            if i>0
                cs(a1,a2) = i;
            else
                % not counted as overlapping entropy
                continue
            end
            

        end
    end
%     keyboard
%     cdsraw = cds;
    cs = Ppair(1).Paa .* cs;
    cds = Ppair(1).Paa .* cds;
elseif NA==3
    for a1=1:Am(1)
        for a2=1:Am(2)
            for a3=1:Am(3)
                s1 = -log2( PA(1).Pa(a1) );
                s2 = -log2( PA(2).Pa(a2) );
                s3 = -log2( PA(3).Pa(a3) );
                sj12 = -log2( Ppair(1).Paa(a1,a2) );
                sj13 = -log2( Ppair(2).Paa(a1,a3) );
                sj23 = -log2( Ppair(3).Paa(a2,a3) );
                sj123 = -log2( Ptrip(1).Paaa(a1,a2,a3) );
                
                % local co-information (entropy overlap)
                i = (sj123 + s1 + s2 + s3 - sj12 - sj13 - sj23);
                % local entropy always positive
                % if local coinformation positive then have overlap
                cds(a1,a2,a3) = i;
                if i>0
                    cs(a1,a2,a3) = i;   
                else
                    % not counted as overlapping entropy
                    continue
                end
            end
        end
    end
    cs = Ptrip(1).Paaa .* cs;
end

% cs
locred = nansum(cs(:));
Hcs = locred;


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
