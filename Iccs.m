function Iccs = Iccs_op_me(A, Pjoint)
% calculate redundancy as from pointwise common change in surprise
% use maximum entropy subejct to pairwise predictor-target marginal
% and all predictor marginal constraints
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

% build distributions for single atom
% Pas = 2d, no maxent solution needed
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

% build pairwise maxent joint distribution for 2 atoms
if NA==2
    thsA = [A{1} A{2}];
    Nv = length(thsA);
    Nv1 = length(A{1});
    Nv2 = length(A{2});
    
    % collapse variables we don't need
    sumover = setdiff(vars, thsA);
    Paas = Pjoint;
    for ii=1:length(sumover)
        Paas = sum(Paas, sumover(ii));
    end
    Paas = squeeze(Paas);
    
    % reorder axes to match order of unique variables in this pair of
    % elements
    % order we want
    Aunq = unique(thsA,'stable');
    % order we have
    [Aunqsrt, Aunqsrtidx] = sort(Aunq);
    % invert order
    [~, Aidx] = sort(Aunqsrtidx);
    Paas = permute(Paas, [Aidx length(Aunq)+1]);
    thsA = changem(thsA, 1:length(Aunq), Aunq);
    Aunq = unique(thsA, 'stable');
    
    % copy duplicate variables as required
    uniquevar_i = 1;
    for allvar_i=1:Nv
        if (uniquevar_i>length(Aunq)) || (thsA(allvar_i) ~= Aunq(uniquevar_i))
            % need to insert a duplicate variable
            var_needed = thsA(allvar_i);
            copy_from = find(thsA==var_needed,1);
            Paas = copy_var(Paas, copy_from, allvar_i);
        else
            % axis order is correct
            uniquevar_i = uniquevar_i + 1;
        end
    end
    
    % joint distribution over all variables
    % in both pairs of elements
    % now should have correct variable axis in correct order
    % collapse A1
    s = size(Paas);
    Paas = reshape(Paas, [prod(s(1:Nv1)) s(Nv1+1:end)]);
    % collapse A2
    s = size(Paas);
    Paas = reshape(Paas, [s(1) prod(s(2:end-1)) s(end)]);
    % pairwise maxent
    P2 = marg_maxent2(Paas);
    Ppair(1).Paas = P2;
    Ppair(1).Paa = squeeze(sum(P2,3));
end

% build triplewise joint element distributions
if NA==3
    thsA = [A{1} A{2} A{3}];
    Nv = length(thsA);
    Nv1 = length(A{1});
    Nv2 = length(A{2});
    Nv3 = length(A{3});

    % collapse variables we don't need
    sumover = setdiff(vars, thsA);
    Paaas = Pjoint;
    for ii=1:length(sumover)
        Paaas = sum(Paaas, sumover(ii));
    end
    Paaas = squeeze(Paaas);

    % reorder axes to match order of unique variables in this pair of
    % elements
    % order we want
    Aunq = unique(thsA,'stable');
    % order we have
    [Aunqsrt, Aunqsrtidx] = sort(Aunq);
    % invert order
    [~, Aidx] = sort(Aunqsrtidx);

    Paaas = permute(Paaas, [Aidx length(Aunq)+1]);
    thsA = changem(thsA, 1:length(Aunq), Aunq);
    Aunq = unique(thsA, 'stable');

    % copy duplicate variables as required
    uniquevar_i = 1;
    for allvar_i=1:Nv
        if (uniquevar_i>length(Aunq)) || (thsA(allvar_i) ~= Aunq(uniquevar_i))
            % need to insert a duplicate variable
            var_needed = thsA(allvar_i);
            copy_from = find(thsA==var_needed,1);
            Paaas = copy_var(Paaas, copy_from, allvar_i);
        else
            % axis order is correct
            uniquevar_i = uniquevar_i + 1;
        end
    end

    % joint distribution over all variables
    % now should have correct variable axes in correct order
    % collapse A1
    s = size(Paaas);
    Nv1 = length(A{1});
    Paaas = reshape(Paaas, [prod(s(1:Nv1)) s(Nv1+1:end)]);
    % collapse A2
    s = size(Paaas);
    Nv2 = length(A{2});
    Paaas = reshape(Paaas, [s(1) prod(s(2:Nv2+1)) s(Nv2+2:end)]);
    % collapse A3
    s = size(Paaas);
    Paaas = reshape(Paaas, [s(1:2) prod(s(3:end-1)) s(end)]);
    
    Pme = marg_maxent_3pred(Paaas);
    Paaas = Pme;
    Ptrip(1).Paaas = Paaas;
    Ptrip(1).Paaa = squeeze(sum(Paaas,4));
    
    % now build pairwise distributions from this maxent solution
    pairs = nchoosek(1:3,2);
    Npair = size(pairs,1);
    Ppair(Npair).Paa = []; % intialize struct
    for pi=1:Npair
        keepax = [pairs(pi,1) pairs(pi,2)];
        % collapse variables we don't need
        sumover = setdiff(1:3, keepax);
        Paas = Paaas;
        for ii=1:length(sumover)
            Paas = sum(Paas, sumover(ii));
        end
        Paas = squeeze(Paas);
        Ppair(pi).Paas = Paas;
        Ppair(pi).Paa = squeeze(sum(Paas,3));
    end
end


% pointwise interaction information
tmp = zeros([Am Sm]);
cds = zeros([Am Sm]);
if NA==1
    for a1=1:Am(1)
        for si=1:Sm
            ds1 = log2( Pele(1).Pas(a1,si) ./ (Pele(1).Pa(a1)*Ps(si)) );
            cds(a1,si) = ds1;
        end
    end
    cds = Pele(1).Pas .* cds;
elseif NA==2
    for a1=1:Am(1)
        for a2=1:Am(2)
            for si=1:Sm
                dsj = (log2( Ppair(1).Paas(a1,a2,si) / (Ppair(1).Paa(a1,a2)*Ps(si)) ));
                ds1 = (log2( Pele(1).Pas(a1,si) ./ (Pele(1).Pa(a1)*Ps(si)) ));
                ds2 = (log2( Pele(2).Pas(a2,si) ./ (Pele(2).Pa(a2)*Ps(si)) ));
                
                num = Pele(1).Pa(a1) * Pele(2).Pa(a2) * Ps(si) * Ppair(1).Paas(a1,a2,si);
                den = Pele(1).Pas(a1,si) * Pele(2).Pas(a2,si) * Ppair(1).Paa(a1,a2);
                ii12 = log2(num ./ den);
                
%                 fprintf(1,'[%d %d %d] : ds1:  %6.3f  ds2:  %6.3f   dsj:  %6.3f   coi:  %6.3f\n',a1,a2,si,ds1,ds2,dsj,-ii12)

                
                overlap = ds1 + ds2 - dsj;
                
                if sign(ds1)==sign(ds2)
                    % change of surprise has same size so possibility of
                    % overlap 
                    if sign(dsj)~=sign(ds1)
%                         fprintf(1,'Warning [%d %d %d] : DSJ sign flip  dsj:  %6.3f  ds1:  %6.3f  ds2:  %6.3f\n',a1,a2,si,dsj,ds1,ds2)
%                         keyboard
                        continue
                    end
                    
                    if sign(overlap)==sign(ds1)
                        % redundant (mis)information
                        if isfinite(overlap) && abs(overlap) > max(abs([ds1 ds2]))
                            fprintf(1,'Warning [%d %d %d] : Overlap larger than individuals. overlap: %6.3f ds1: %6.3f  ds2: %6.3f\n',a1,a2,si,overlap,ds1,ds2)
                        end
                        cds(a1,a2,si) = overlap;
                    end
                end   
            end
        end
    end
    cds = Ppair(1).Paas .* cds;
elseif NA==3
    for a1=1:Am(1)
        for a2=1:Am(2)
            for a3=1:Am(3)
                for si=1:Sm
                    ds123 = log2( Ptrip(1).Paaas(a1,a2,a3,si) / (Ptrip(1).Paaa(a1,a2,a3)*Ps(si)) );
                    ds1 = log2( Pele(1).Pas(a1,si) ./ (Pele(1).Pa(a1)*Ps(si)) );
                    ds2 = log2( Pele(2).Pas(a2,si) ./ (Pele(2).Pa(a2)*Ps(si)) );
                    ds3 = log2( Pele(3).Pas(a3,si) ./ (Pele(3).Pa(a3)*Ps(si)) );
                    ds12 = (log2( Ppair(1).Paas(a1,a2,si) / (Ppair(1).Paa(a1,a2)*Ps(si)) ));
                    ds13 = (log2( Ppair(2).Paas(a1,a3,si) / (Ppair(2).Paa(a1,a3)*Ps(si)) ));
                    ds23 = (log2( Ppair(3).Paas(a2,a3,si) / (Ppair(3).Paa(a2,a3)*Ps(si)) ));

                    if (sign(ds1)==sign(ds2)) && (sign(ds2)==sign(ds3)) 
                        % change of surprise has same sign for all 3
                        % variables, so possibility of overlap
                        if sign(ds123)~=sign(ds1)
                            
%                         if sign(ds123)~=sign(ds1) || sign(ds12)~=sign(ds1) || sign(ds13) ~=sign(ds1) || sign(ds23) ~= sign(ds1)
%                             fprintf(1,'Warning [%d %d %d %d] : DSJ sign flip  dsj: %6.3f  ds1: %6.3f  ds2: %6.3f ds3: %6.3f\n',a1,a2,a3,si,ds123,ds1,ds2,ds3)
                            continue
                        end
                        overlap = ds1 + ds2 + ds3 - ds12 - ds13 - ds23 + ds123;
                        if sign(overlap)==sign(ds1)
                            % redundant (mis)information
                            if isfinite(overlap) && abs(overlap) > max(abs([ds1 ds2 ds3]))
                                fprintf(1,'Warning [%d %d %d %d] : Overlap larger than individual overlap: %6.3f ds1: %6.3f  ds2: %6.3f  ds3: %6.3f\n',a1,a2,a3,si,overlap,ds1,ds2,ds3)
                            end
                            cds(a1,a2,a3,si) = overlap;
                        end
                    end
%                     tmp(a1,a2,a3,si) = ds123;

                end
            end
        end
    end
    cds = Ptrip(1).Paaas .* cds;
end
% cds
locred = nansum(cds(:));
Iccs = locred;

function y = fixsign(x)
% fix the sign of things close to zero following dit (which uses
% np.isclose)
% absolute(a - b) <= (atol + rtol * absolute(b))
atol = 1e-8;
if abs(x) <= atol
    y = 0.0;
else
    y = x;
end

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



