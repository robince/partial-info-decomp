function Iccs = Iccs_Pind(A, Pjoint)
% calculate redundancy between a set of sources
% from pointwise common change in surprise
% using conditional independent joint-source distribution
%
% A - cell array of sources
% Pjoint - full joint distribution (S last axis)

isclosefp = @(a,b) abs(a - b) <= eps(max(abs(a), abs(b)));
s = size(Pjoint);
Sm = s(end); % number of target values
Nx = length(s)-1; % number of dependent variables
vars = 1:Nx;

NA = length(A);
if NA>3
    error('Iccs: only 3 sources supported')
end
PA(NA).Pa = []; % intialize struct
Am = zeros(1,NA); % number of symbols in each element

% Ps
Ps = Pjoint;
for xi=1:Nx
    Ps = squeeze(sum(Ps,1));
end

% sort variables within each source
A = cellfun(@sort, A, 'Unif',false);

% build distributions for each source
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
    PA(ai).Pas = Pas;
    % unconditional distribution P(a)
    PA(ai).Pa = squeeze(sum(Pas,2));
    % conditional distribution P(a|s)
    PA(ai).Pacs = bsxfun(@rdivide, PA(ai).Pas, Ps);
    Am(ai) = size(Pas,1);
end

% build pairwise joint source distributions
if NA==2
    pairs = nchoosek(1:NA,2);
    Npair = size(pairs,1);
    Ppair(Npair).Paa = []; % intialize struct
    for pi=1:Npair
        p1 = pairs(pi,1);
        p2 = pairs(pi,2);
        Ppair(pi).Paacs = zeros(Am(p1),Am(p2),Sm);
        for si=1:Sm
            Paacs(:,:,si) = PA(p1).Pacs(:,si) * PA(p2).Pacs(:,si)';
        end
        Ppair(pi).Paacs = Paacs;
        Paas = bsxfun(@times, Paacs, reshape(Ps,[1 1 Sm]));
        Ppair(pi).Paas = Paas;
        Ppair(pi).Paa = nansum(Paas,3);
    end
end

% build triplewise joint source distributions
if NA==3
    Paaacs = zeros(Am(1),Am(2),Am(3),Sm);
    for si=1:Sm
        Paac = PA(1).Pacs(:,si) * PA(2).Pacs(:,si)';
        Paaacs(:,:,:,si) = bsxfun(@times, Paac, reshape(PA(3).Pacs(:,si), [1 1 Am(3)]));
    end
    Ptrip(1).Paaacs = Paaacs;
    Paaas = bsxfun(@times, Paaacs, reshape(Ps,[1 1 1 Sm]));
    Ptrip(1).Paaas = Paaas;
    Ptrip(1).Paaa = nansum(Paaas,4);
    
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
cds = zeros([Am Sm]);
if NA==1
    for a1=1:Am(1)
        for si=1:Sm
            ds1 = log2( PA(1).Pas(a1,si) ./ (PA(1).Pa(a1)*Ps(si)) );
            cds(a1,si) = ds1;
        end
    end
%     keyboard
    cds = PA(1).Pas .* cds;
elseif NA==2
    for a1=1:Am(1)
        for a2=1:Am(2)
            for si=1:Sm
                dsj = log2( Ppair(1).Paas(a1,a2,si) / (Ppair(1).Paa(a1,a2)*Ps(si)) );
                ds1 = log2( PA(1).Pas(a1,si) ./ (PA(1).Pa(a1)*Ps(si)) );
                ds2 = log2( PA(2).Pas(a2,si) ./ (PA(2).Pa(a2)*Ps(si)) );

                overlap = ds1 + ds2 - dsj;
                if isfinite(overlap) && (abs(overlap)-min(abs([ds1 ds2])))>2*eps(abs(overlap))
                    fprintf(1,'Warning [%d %d %d] : Overlap larger than individuals. overlap: %6.3f ds1: %6.3f  ds2: %6.3f\n',a1,a2,si,overlap,ds1,ds2)
                    fprintf(1,'overlap: %d ds1: %d  ds2: %d\n',sign(overlap),sign(ds1),sign(ds2))
                    
                end
                if sign(ds1)==sign(ds2)
                    % change of surprise has same size so possibility of
                    % overlap
                    if sign(dsj)~=sign(ds1)
%                         fprintf(1,'Warning [%d %d %d] : DSJ sign flip  dsj:  %6.3f  ds1:  %6.3f  ds2:  %6.3f\n',a1,a2,si,dsj,ds1,ds2)
%                         keyboard
                        continue
                    end
                    overlap = ds1 + ds2 - dsj;
                    
                    if sign(overlap)==sign(ds1)
                        % redundant (mis)information
                        if isfinite(overlap) && (abs(overlap)-min(abs([ds1 ds2])))>2*eps(abs(overlap))
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
                    ds1 = log2( PA(1).Pas(a1,si) ./ (PA(1).Pa(a1)*Ps(si)) );
                    ds2 = log2( PA(2).Pas(a2,si) ./ (PA(2).Pa(a2)*Ps(si)) );
                    ds3 = log2( PA(3).Pas(a3,si) ./ (PA(3).Pa(a3)*Ps(si)) );
                    ds12 = log2( Ppair(1).Paas(a1,a2,si) / (Ppair(1).Paa(a1,a2)*Ps(si)) );
                    ds13 = log2( Ppair(2).Paas(a1,a3,si) / (Ppair(2).Paa(a1,a3)*Ps(si)) );
                    ds23 = log2( Ppair(3).Paas(a2,a3,si) / (Ppair(3).Paa(a2,a3)*Ps(si)) );
                    
                    if (sign(ds1)==sign(ds2)) && (sign(ds2)==sign(ds3))
                        % change of surprise has same sign for all 3
                        % variables, so possibility of overlap
                        if sign(ds123)~=sign(ds1)
%                             fprintf(1,'Warning [%d %d %d %d] : DSJ sign flip  dsj: %6.3f  ds1: %6.3f  ds2: %6.3f ds3: %6.3f\n',a1,a2,a3,si,ds123,ds1,ds2,ds3)
                            continue
                        end
                        overlap = ds1 + ds2 + ds3 - ds12 - ds13 - ds23 + ds123;
                        if sign(overlap)==sign(ds1)
                            % redundant (mis)information                            
                            if isfinite(overlap) && (abs(overlap)-max(abs([ds1 ds2 ds3])))>6*eps(min(abs([overlap ds1 ds2 ds3])))
                                fprintf(1,'Warning [%d %d %d %d] : Overlap larger than individual overlap: %6.3f ds1: %6.3f  ds2: %6.3f  ds3: %6.3f\n',a1,a2,a3,si,overlap,ds1,ds2,ds3)
                            end
                            cds(a1,a2,a3,si) = overlap;
                        end
                    end
                end
            end
        end
    end
    cds = Ptrip(1).Paaas .* cds;
end

locred = nansum(cds(:));
Iccs = locred;



